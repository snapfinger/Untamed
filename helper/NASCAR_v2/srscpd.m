%
% result = srscpd(TS, R, option)
% 
% Description:
%     srscpd framework (N-D supported)
% 
% Input:
%     TS - tensor
%     R - desired rank
%     option - option structure, see code
% 
% Output:
%     result - decomposition result struct with rank from 1 to R
% 
% Copyright:
%     2017-2020 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     2.5.0
% Date:
%     2021/10/31
%

function result = srscpd(TS, R, option)
    
    if nargin == 0
        result = srscpd('als');
        return;
    elseif nargin == 1
        if ischar(TS) || isa(TS, 'function_handle')
            option = struct;
            option.isStats = false;
            option.maxNumFitRes = 10;
            option.isVerbose = true;
            option.nonnegative = [];
            option.cacheMTS = true;
            option.rankOneOptALS = cpALS();
            option.rankOneInit = 'random';
            if ischar(TS)
                switch TS
                    case 'als'
                        option.alg = 'als';
                        option.optAlg = cpALS();
                        option.algFunc = @cpALS;
                    case 'opt'
                        option.alg = 'opt';
                        option.optAlg = cpOpt();
                        option.algFunc = @cpOpt;
                    otherwise
                        warning('unknown algorithm, use ALS');
                        option.alg = 'als';
                        option.optAlg = cpALS();
                        option.algFunc = @cpALS;
                end
            elseif isa(TS, 'function_handle')
                option.alg = 'custom';
                option.algFunc = TS;
                option.optAlg = TS();
            end
            option.saveToFile = [];
            option.logFile = [];
            option.resumeFrom = [];
            result = option;
            return;
        end
    end
    
    if ~exist('option', 'var') || isempty(option)
        option = srscpd();
    end
    
    isVerbose = option.isVerbose;
    
    % rank 1 ALS option
    rankOneInit = option.rankOneInit;
    optALS = option.rankOneOptALS;
    if iscell(rankOneInit)
        optALS.init = rankOneInit{1};
    elseif ischar(rankOneInit)
        optALS.init = rankOneInit;
    end
    optALS.nonnegative = option.nonnegative;
    optALS.cacheMTS = option.cacheMTS;
    if ~isVerbose, optALS.printItv = 0; end
    
    % rank r option
    algFunc = option.algFunc;
    optAlg = option.optAlg;
    optAlg.nonnegative = option.nonnegative;
    optAlg.cacheMTS = option.cacheMTS;
    if ~isVerbose
        optAlg.optSolver.isVerbose = false;
        optAlg.optSolver.printItv = 0;
    end
    
    % start srscpd framework
    if isempty(option.logFile)
        if ~isempty(option.saveToFile)
            [savePath, saveName, ~] = fileparts(option.saveToFile);
            logFile = fullfile(savePath, [saveName '.log']);
            isLog = true;
        else
            isLog = false;
        end
    else
        logFile = option.logFile;
        isLog = true;
    end
    
    if isLog
        if exist(logFile, 'file') && isempty(option.resumeFrom)
            delete(logFile);
        end
        diary(logFile);
    end
    
    ts = tic;
    if isVerbose
        fprintf('\n========== srscpd started ==========\n');
        fprintf('starting time: %s\n', datestr(now));
    end
    
    % if start from from 1
    if isempty(option.resumeFrom)
        if isVerbose, fprintf('\nfresh start from rank 1\n\n'); end
        
        result  = struct;
        
        if isVerbose, disp('fit rank 1 tensor as basis'); end
        % fit the first rank 1 tensor
        [U, lambda, output] = cpALS(TS, 1, optALS);
        % re-fit if failed
        c = 1;
        while (~output.Flag) && (c < option.maxNumFitRes)
            c = c + 1;
            if isVerbose, disp(['R = 1 failed, try again (' num2str(c) ')']); end
            [U, lambda, output] = cpALS(TS, 1, optALS);
        end
        % if failed too many times at the first round, just quit
        if ~output.Flag
            if isVerbose, disp('still failed to fit the first rank 1 tensor, quit'); end
            return;
        end

        result(1).U = U;
        result(1).Lambda = lambda;
        result(1).Output = output;

        % calculate statistics if needed
        if option.isStats
            result(1).Stats.EV = output.EV;
            result(1).Stats.CCD = CCDWrapper(TS, U, lambda);
        end
        
        r = 1;
    else
        % if resume from certain rank
        result = option.resumeFrom;
        r = length(result);
        U = result(r).U;
        lambda = result(r).Lambda;
        if isVerbose, fprintf('\nresume from rank %d\n', r); end
    end
    
    if ~isempty(option.saveToFile)
        save(option.saveToFile, 'result');
    end
    
    % calculate the residue
    TSRes = TS - cpFull(U, lambda);
    
    % iterative over the rest of ranks
    N = length(U);
    
    for m = r+1:R
        if isVerbose, fprintf('\nfit rank 1 tensor to residue as part of warm start\n'); end
        
        % fit rank 1 tensor to the residue
        if iscell(rankOneInit), optALS.init = rankOneInit{m}; end
        [URes, lambdaRes, output] = cpALS(TSRes, 1, optALS);
        
        % re-fit if failed
        c = 1;
        while (~output.Flag) && (c < option.maxNumFitRes)
            c = c + 1;
            if isVerbose, disp(['did not find good init from residue, try again (' num2str(c) ')']); end
            [URes, lambdaRes, output] = cpALS(TSRes, 1, optALS);
        end
        
        % just use random if still failed
        if ~output.Flag
            if isVerbose, disp('still could not find good init, just use random'); end
            URes = cell(N, 1);
            for n = 1:N
                URes{n} = rand(size(TS, n), 1);
            end
        end
        
        % equally spread the scale lambda to each dimension
        UInit = cell(N, 1);
        for n = 1:N
            a = U{n};
            b = URes{n};
            UInit{n} = [a .* (lambda.^(1/N))', b .* (lambdaRes.^(1/N))'];
        end
        
        % fit rank r tensor
        if isVerbose, fprintf('\nfit rank %d tensor\n', m); end
        
        optAlg.init = UInit;
        [U, lambda, output] = algFunc(TS, m, optAlg);
        
        result(m).U = U;
        result(m).Lambda = lambda;
        result(m).Output = output;
        
        % calculate statistics if needed
        if option.isStats
            result(m).Stats.EV = output.EV;
            result(m).Stats.CCD = CCDWrapper(TS, U, lambda);
        end
        
        % calculate the residue
        TSRes = TS - cpFull(U, lambda);
        
        if ~isempty(option.saveToFile)
            save(option.saveToFile, 'result');
        end
    end
    
    te = toc(ts);
    if isVerbose
        fprintf('\nfinish time: %s\n', datestr(now));
        te2 = te ./ (60*60*24);
        fprintf('time elapsed: %s day %s hr %s min %s sec\n', ...
                datestr(te2, 'dd'), datestr(te2, 'HH'), ...
                datestr(te2, 'MM'), datestr(te2, 'SS'));
        fprintf('\n========== srscpd finished ==========\n');
    end
    
    if isLog
        diary('off');
    end
end
