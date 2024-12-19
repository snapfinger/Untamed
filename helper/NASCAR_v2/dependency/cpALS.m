%
% [U, lambda, output] = cpALS(TS, R, option)
% 
% Description:
%     Alternating Least Square (ALS) algorithm for CP decomposition (N-D supported)
% 
% Input:
%     TS - N-way tensor
%     R - desired rank
%     option - option structure, see code
% 
% Output:
%     U - the components
%     lambda - scale corrsponding to the components
%     output - auxiliary output
% 
% Copyright:
%     2017-2020 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     2.4.3
% Date:
%     2021/07/04
%

function [U, lambda, output] = cpALS(TS, R, option)
    
    if nargin == 0
        option = struct;
        
        % initialzation, input is either a matrix or 'random' or 'svd'
        option.init = 'random';
        
        % the first iterating dimensions
        option.firstItrDims = [];
        
        % the overall iterating dimensions
        % (fixed dimensions should be the complement of this set from the whole set 1:N)
        option.itrDims = []; 
        
        % convergence criterion
        option.tol = 1e-5;
        % max number of iterations
        option.maxNumItr = 100;
        
        % whether to cache the matricized tensor or not
        % caching will speed up the computation but require more memory
        option.cacheMTS = true;
        
        % constrained dimensions
        % 0 - no constraint, 1 - smoothness, 2 - sparsity 3 - TV
        option.const = [];
        % huber approx eps for l1 norm (sparsity)
        option.epsHuber = 1e-3;
        % regurlarization parameters for each dimension
        option.regParam = [];
        % non-negative constraint
        option.nonnegative = [];
        
        % output the history of cost function value
        option.isCostFV = false;
        % output the history of the convergence
        option.isFC = false;
        
        % print every N iterations (first 5 iterations are always printed)
        option.printItv = 10;
        
        % options for TFOCS
        if exist('tfocs', 'file')
            option.useTFOCS = true;
            optTFOCS = tfocs;
            optTFOCS.tol = 1e-6;
            optTFOCS.restart = 5;
            optTFOCS.maxIts = 100;
            optTFOCS.alg = 'AT';
            optTFOCS.printEvery = 0;
            optTFOCS.debug = false;
            option.optionTFOCS = optTFOCS;
        else
            option.useTFOCS = false;
            option.optionTFOCS = [];
        end
        
        U = option;
        lambda = [];
        output = [];
        return;
    end
    
    warning('off', 'all');
    
    % initialize options
    % need to define dimensions for initializations
    N = ndims(TS);
    sz = size(TS);
    clsName = class(TS);
    
    if option.useTFOCS
        if exist('tfocs', 'file') && isa(TS, 'double')
            isUseTFOCS = true;
        else
            if ~exist('tfocs', 'file'), disp('Missing tfocs toolbox'); end
            if ~isa(TS, 'double'), disp('Need double precision to use the sparse matrix in tfocs'); end
            disp('Use mldivide instead')
            isUseTFOCS = false;
        end
    else
        isUseTFOCS = false;
    end
    
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
    
    if isempty(option.firstItrDims)
        firstItrDims = 1:N;
    else
        firstItrDims = option.firstItrDims;
    end
    
    if isempty(option.itrDims)
        itrDims = 1:N;
    else
        itrDims = option.itrDims;
    end
    
    if isempty(option.const)
        const = zeros(N, 1);
    else
        if length(option.const) ~= N
            error('regularization type dimension mismatches the data');
        end
        const = option.const;
    end
    
    if (~isUseTFOCS) && any(const)
        disp('Regularization is not yet supported without TFOCS');
        U = [];
        lambda = [];
        output = [];
        return;
    end
   
    if isempty(option.regParam)
        mu = zeros(N, 1);
    else
        if length(option.regParam) ~= N
            error('regularization parameter dimension mismatches the data');
        end
        mu = option.regParam;
    end
    
    if isempty(option.nonnegative)
        nonneg = false(N, 1);
    else
        if length(option.nonnegative) ~= N
            error('nonnegativity constraint dimension mismatches the data');
        end
        nonneg = option.nonnegative;
    end
    
    maxNumItr = option.maxNumItr;
    epsHuber = option.epsHuber;
    isCostFV = option.isCostFV;
    isFC = option.isFC;
    printItv = option.printItv;
    
    % ALS algorithm starts
    U = UInit;
    
    if isCostFV, CostFV = zeros(maxNumItr, 1, clsName); end
    if isFC, FC = zeros(maxNumItr, 1, clsName); end
    
    V = zeros(R, R, N, clsName);
    for m = 1:N
        V(:, :, m) = (U{m}') * U{m};
    end
    
    for m = 1:maxNumItr
        UOld = U;
        
        for n = firstItrDims % iterate over modes
            idx = [1:n-1 n+1:N];
            
            A = (prod(V(:, :, idx), 3));
            if option.cacheMTS
                B = (MTS{n} * krProd(U(idx(end:-1:1))))';
            else
                B = (matricize(TS, n) * krProd(U(idx(end:-1:1))))';
            end
            
            if isUseTFOCS
                A2 = nBlockDiag(A, sz(n));
                B2 = B(:);
                x0 = U{n}';
                x0 = x0(:);

                if const(n) == 0
                    smoothF = smooth_quad;
                    linearF = {A2, -B2};
                elseif const(n) == 1 % smoothness
                    smoothF2 = tfunc_scale(smooth_quad, mu(n), 1, 0);
                    smoothF = {smooth_quad, smoothF2};
                    D = linop_diff_h([R, sz(n)]);
                    linearF = {A2, -B2; D, 0};
                elseif const(n) == 2 % sparsity
                    smoothF2 = tfunc_scale(smooth_huber(epsHuber), mu(n), 1, 0);
                    smoothF = {smooth_quad, smoothF2};
                    linearF = {A2, -B2; 1, 0};
                elseif const(n) == 3 % total variantion
                    smoothF2 = tfunc_scale(smooth_huber(epsHuber), mu(n), 1, 0);
                    smoothF = {smooth_quad, smoothF2};
                    D = linop_diff_h([R, sz(n)]);
                    linearF = {A2, -B2; D, 0};
                end

                if nonneg(n)
                    [X2, ~] = tfocs(smoothF, linearF, proj_Rplus, x0, option.optionTFOCS);
                else
                    [X2, ~] = tfocs(smoothF, linearF, [], x0, option.optionTFOCS);
                end
                X = reshape(X2, R, sz(n));
            else
                X = A \ B;
                if nonneg(n)
                    X(X < 0) = 0;
                end
            end
            
            lambda = sqrt(diag(X * X'));
            X = bsxfun(@rdivide, X, lambda);

            U{n} = X';
            V(:, :, n) = X * X';
        end
        
        % check convergence (relative fro norm of diff)
        [~, xDiff] = cpDiffX(U, UOld, 2);
        if isFC, FC(m) = xDiff; end
        
        % print
        if (printItv > 0) && ((m <= 5) || (mod(m, printItv) == 0))
            fprintf('%d: dx = %.3e\n', m, xDiff);
        end
        
        % calculate cost function value if it is requested
        if isCostFV
            [dfFro, ~, ~] = cpDiff(TS, U, lambda);
            CostFV(m) = dfFro;
        end
        
        % reset to full or constrained modes after first itr
        firstItrDims = itrDims;
        
        if (m > 2) && (xDiff < option.tol)
            break;
        end
    end
    
    clearvars MTS;
    
    [dfFro, ~, EV] = cpDiff(TS, U, lambda);

    if isnan(dfFro)
        if printItv > 0
            disp('invalid result');
        end
        output = struct;
        output.Flag = false;
        return;
    end
    
    if printItv > 0
        if m == maxNumItr
            disp('reached the max number of iterations');
            fprintf('df = %.2f, dx = %.2e, EV = %.4f\n', dfFro, xDiff, EV);
        else
            fprintf('Converged: df = %.2f, dx = %.2e, EV = %.4f\n', dfFro, xDiff, EV);
        end
    end
    
    output = struct;
    output.Flag = true;
    output.numItr = m;
    output.EV = EV;
    
    if isCostFV, output.CostFV = CostFV(1:m); end
    if isFC, output.FC = FC(1:m); end
end
