%
% [U, lambda, output] = cpOpt(TS, R, option)
% 
% Description:
%     optimization-based CP decomposition (N-D supported)
% 
% Input:
%     TS - tensor
%     R - desired rank
%     option - option structure, see code
% 
% Output:
%     U  - the components
%     lambda - scale corrsponding to the components
%     output - auxiliary output
% 
% Copyright:
%     2017-2021 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.2.1
% Date:
%     2021/10/31
%

function [U, lambda, output] = cpOpt(TS, R, option)

    if nargin == 0
        option = struct;
        
        % initialzation
        option.init = 'random';
        option.normRegParam = 0.01;
        
        option.cacheMTS = true;
        
        % non-negative constraint
        option.nonnegative = [];
        
        option.solver = @adam_solver;
        optSolver = option.solver();
        option.optSolver = optSolver;
        
        U = option;
        lambda = [];
        output = [];
        return;
    end
    
    if ~exist('option', 'var') || isempty(option)
        option = cpOpt();
    end
    
    N = ndims(TS);
    sz = size(TS);
    clsName = class(TS);
    
    % set initialization
    if ischar(option.init)
        UInit = cell(N, 1);
        for m = 1:N
            if strcmpi(option.init, 'random')
                UInit{m} = rand(size(TS, m), R, clsName);
            elseif strcmpi(option.init, 'svd')
                [LSV, ~, ~] = svds(matricize(double(TS), m), R);
                UInit{m} = cast(LSV, clsName);
            else
                error('initialization not specified');
            end
        end
    else
        UInit = option.init;
        if ~isa(UInit{1}, clsName)
            for m = 1:N
                UInit{m} = cast(UInit{m}, clsName);
            end
        end
    end
    
    if option.cacheMTS
        MTS = cell(N, 1);
        for m = 1:N
            MTS{m} = matricize(TS, m);
        end
    else
        MTS = {};
    end
    
    if isempty(option.nonnegative)
        nonneg = false(N, 1);
    else
        if length(option.nonnegative) ~= N
            error('nonnegativity constraint dimension mismatches the data');
        end
        nonneg = option.nonnegative;
    end
    
    nonnegMask = cell(N, 1);
    for m = 1:N
        if nonneg(m)
            nonnegMask{m} = true(size(UInit{m}));
        else
            nonnegMask{m} = false(size(UInit{m}));
        end
    end
    nonnegMask = U2Uvec(nonnegMask);
    option.optSolver.nonnegMask = nonnegMask;
    
    normRegLmd = option.normRegParam;
    
    objFun = @cpCostFunction;
    UInitVec = U2Uvec(UInit);
    
    [UVec, output] = option.solver(objFun, UInitVec, option.optSolver);
    U = Uvec2U(UVec);
    
    clearvars MTS;
    
    % store norm in lambda
    lambda = ones(R, 1, clsName);
    for m = 1:N
        t = sqrt(diag(U{m}' * U{m}));
        lambda = lambda .* t;
        U{m} = bsxfun(@rdivide, U{m}', t)';
    end
    
    [dfFro, ~, EV] = cpDiff(TS, U, lambda);
    
    if isnan(dfFro)
        output.Flag = false;
    else
        output.EV = EV;
        output.Flag = true;
    end
    
%     fprintf('EV = %.4f\n', EV);
    
    function [funcValue, UGradvec] = cpCostFunction(Uvec2, isFVal)
        if ~exist('isFVal', 'var') || isempty(isFVal)
            isFVal = true;
        end
        
        U2 = Uvec2U(Uvec2);

        % calculate function value
        if isFVal
            if isempty(MTS)
                B = matricize(TS, 1) - U2{1} * krProd(U2(end:-1:2))';
            else
                B = MTS{1} - U2{1} * krProd(U2(end:-1:2))';
            end
            UNorm2 = zeros(1, clsName);
            for n = 1:N
                UNorm2 = UNorm2 + sum(diag(U2{n}' * U2{n}));
            end
            funcValue = 0.5 * (norm(B, 'fro')^2 + normRegLmd * UNorm2);
%             f1 = norm(B, 'fro')^2
%             f2 = UNorm2
        else
            funcValue = [];
        end

        % calculate gradient
        UGrad = cell(size(U2));

        V = zeros(R, R, N, clsName);
        for n = 1:N
            V(:, :, n) = (U2{n}') * U2{n};
        end

        for n = 1:N
            idx = [1:n-1 n+1:N];
            A = (prod(V(:, :, idx), 3));
            if isempty(MTS)
                B = matricize(TS, n) * krProd(U2(idx(end:-1:1)));
            else
                B = MTS{n} * krProd(U2(idx(end:-1:1)));
            end
            UGrad{n} = -B + U2{n} * A + normRegLmd * U2{n};
        end

        UGradvec = U2Uvec(UGrad);
    end

    function Uvec = U2Uvec(U)
        blk = sz * R;
        blkStep = [0, cumsum(blk)]';

        Uvec = zeros(sum(blk), 1, clsName);
        for n = 1:N
            Uvec(blkStep(n)+1:blkStep(n+1)) = reshape(U{n}, blk(n), 1);
        end
    end

    function U = Uvec2U(Uvec)
        blk = sz * R;
        blkStep = [0, cumsum(blk)]';

        U = cell(N, 1);
        for n = 1:N
            U{n} = reshape(Uvec(blkStep(n)+1:blkStep(n+1)), sz(n), R);
        end
    end

end

