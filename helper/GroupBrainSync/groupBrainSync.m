%
% [Y, O] = groupBrainSync(data, option)
% 
% Description:
%     group BrainSync algorithm
% 
% Input:
%     data - T x V x S data tensor
%     option - option structure
% 
% Output:
%     Y - group atlas
%     O - orthogonal transformation matrices
% 
% Copyright:
%     2019-2021 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian Li (Andrew), Haleh Akrami, Anand A. Joshi
% Revision:
%     2.1.1
% Date:
%     2020/03/30
%

function [Y, O] = groupBrainSync(data, option)

    if nargin == 0
        option = struct;
        option.isDataMatrix = true;
        option.dataLoadingFunc = [];
        option.batchSize = 1;
        option.numCPU = 0;
        option.epsX = 1e-3;
        option.isVerbose = true;
        option.save_iter = 20;
        Y = option;
        return;
    end
    
    if ~exist('option', 'var') || isempty(option)
        option = groupBrainSync();
    end
    
    isMat = option.isDataMatrix;
    DLFunc = option.dataLoadingFunc;
    bSize = option.batchSize;
    
    isPar = false;
    if option.numCPU > 0
        if ~parpoolOperator('isopen')
            parpoolOperator('open', option.numCPU);
        end
        isPar = true;
    end    
    isVerbose = option.isVerbose;
    
    if isMat
        [T, V, S] = size(data);
        clsName = class(data);
        X = zeros(T, V, S, clsName);
    else
        data2 = DLFunc(data{1});
        [T, V] = size(data2);
        clsName = class(data2);
        S = length(data);
        numDL = ceil(S / bSize);
    end
    
    if isVerbose, disp('initialize Os and Y'); end
    O = zeros(T, T, S, clsName);
    Y = zeros(T, V, clsName);
    
    if isMat
        if isPar
            if isVerbose, parProgressTracker('s', S); end
            parfor m = 1:S
                R = 2 * rand(T, T) - 1;
                R = (R * R')^(-1/2) * R;
                O(:, :, m) = R';
                X(:, :, m) = O(:, :, m) * data(:, :, m);
                if isVerbose, parProgressTracker('p'); end
            end
            if isVerbose, parProgressTracker('e'); end
        else
            if isVerbose, strLen = progressTracker(0, S, 0, 50); end
            for m = 1:S
                R = 2 * rand(T, T) - 1;
                R = (R * R')^(-1/2) * R;
                O(:, :, m) = R';
                X(:, :, m) = O(:, :, m) * data(:, :, m);
                if isVerbose, strLen = progressTracker(m, S, strLen, 50); end
            end
            clear R;
        end
        Y = mean(X, 3);
    else
        if isVerbose, strLen = progressTracker(0, numDL, 0, 50); end
        for n = 1:numDL
            X = zeros(T, V, bSize, clsName);
            if isPar
                O2 = zeros(T, T, bSize, clsName);
                parfor m = 1:bSize
                    k = (n - 1) * bSize + m;
                    if k <= S
                        R = 2 * rand(T, T) - 1;
                        R = (R * R')^(-1/2) * R;
                        O2(:, :, m) = R';
                        data2 = DLFunc(data{k});
                        X(:, :, m) = O2(:, :, m) * data2;
                    end
                end
                idxS = (n - 1) * bSize + 1;
                idxE = min(n * bSize, S);
                idxSE2 = (1:bSize) + (n - 1) * bSize <= S;
                O(:, :, idxS:idxE) = O2(:, :, idxSE2);
                clear O2;
            else
                for m = 1:bSize
                    k = (n - 1) * bSize + m;
                    if k <= S
                        R = 2 * rand(T, T) - 1;
                        R = (R * R')^(-1/2) * R;
                        O(:, :, k) = R';
                        data2 = DLFunc(data{k});
                        X(:, :, m) = O(:, :, k) * data2;
                    end
                end
                clear R data2;
            end
            Y = Y + sum(X, 3);
            if isVerbose, strLen = progressTracker(n, numDL, strLen, 50); end
        end
        Y = Y / S;
    end
    
    err = 1e10; YOld = Y;
    c = 0;
    
    if isVerbose, disp('start iterations'); end
    while err > option.epsX
        if isMat
            if isPar
                if isVerbose, parProgressTracker('s', S); end
                parfor m = 1:S
                    C = Y * data(:, :, m)';
                    [U2, ~, V2] = svd(C);
                    O(:, :, m) = U2 * V2';
                    X(:, :, m) = O(:, :, m) * data(:, :, m);
                    if isVerbose, parProgressTracker('p'); end
                end
                if isVerbose, parProgressTracker('e'); end
            else
                if isVerbose, strLen = progressTracker(0, S, 0, 50); end
                for m = 1:S
                    C = Y * data(:, :, m)';
                    [U2, ~, V2] = svd(C);
                    O(:, :, m) = U2 * V2';
                    X(:, :, m) = O(:, :, m) * data(:, :, m);
                    if isVerbose, strLen = progressTracker(m, S, strLen, 50); end
                end
                clear C U2 V2;
            end
            Y = mean(X, 3);
        else
            if isVerbose, strLen = progressTracker(0, numDL, 0, 50); end
            Y2 = zeros(T, V, clsName);
            for n = 1:numDL
                X = zeros(T, V, bSize, clsName);
                if isPar
                    O2 = zeros(T, T, bSize, clsName);
                    parfor m = 1:bSize
                        k = (n - 1) * bSize + m;
                        if k <= S
                            data2 = DLFunc(data{k});
                            C = Y * data2';
                            [U2, ~, V2] = svd(C);
                            O2(:, :, m) = U2 * V2';
                            X(:, :, m) = O2(:, :, m) * data2;
                        end
                    end
                    idxS = (n - 1) * bSize + 1;
                    idxE = min(n * bSize, S);
                    idxSE2 = (1:bSize) + (n - 1) * bSize <= S;
                    O(:, :, idxS:idxE) = O2(:, :, idxSE2);
                    clear O2;
                else
                    for m = 1:bSize
                        k = (n - 1) * bSize + m;
                        if k <= S
                            data2 = DLFunc(data{k});
                            C = Y * data2';
                            [U2, ~, V2] = svd(C);
                            O(:, :, k) = U2 * V2';
                            X(:, :, m) = O(:, :, k) * data2;
                        end
                    end
                    clear C U2 V2 data2;
                end
                Y2 = Y2 + sum(X, 3);
                if isVerbose, strLen = progressTracker(n, numDL, strLen, 50); end
            end
            Y = Y2 / S;
        end
        
        YDiff = Y - YOld;
        err = norm(YDiff(:)) / norm(YOld(:));
        YOld = Y;
        
        if isVerbose, fprintf('\niter %d: rel err: %.4g\n\n', c, err); end
        
        if ~mod(c, option.save_iter)
            save(['GSync_Atlas1_' num2str(c) '_err' num2str(err) '.mat' ], 'Y');
        end
        
        c = c + 1;
        
        
        
    end
end