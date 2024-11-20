%
% [x, output] = adam_solver(objFun, x0, option)
% 
% Description:
%     Epoch-based adaptive moment estimation solver (Adam) with optional Nesterov acceleration
% 
% Input:
%     objFun - objective function with output function value and gradient
%     x0 - starting point
%     option - option structure, see code
% 
% Output:
%     x  - solution at convergence
%     output - auxiliary output
% 
% Copyright:
%     2017-2021 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.1.6
% Date:
%     2021/10/31
%

function [x, output] = adam_solver(objFun, x0, option)

    if nargin == 0
        option = struct;
        option.learningRate = 0.001;
        option.beta1 = 0.9;
        option.beta2 = 0.999;
        option.maxNumIter = 1e4;
        option.tolX = 1e-4;
        option.tolFunc = 1e-6;
        option.epsilon = sqrt(eps);
        option.isNesterovAccelerated = true;
        option.printItv = 10;
        option.nonnegMask = [];
        option.epoch = 10;
        option.isVerbose = true;
        
        x = option;
        output = [];
        return;
    end
    
    if (~exist('option', 'var') || isempty(option))
        option = adam_solver();
    end
    
    clsName = class(x0);
    
    maxNumIter = option.maxNumIter;
    eta = option.learningRate;
    b1 = option.beta1;
    b2 = option.beta2;
    epsilon = option.epsilon;
    printItv = option.printItv;
    epoch = option.epoch;
    isVerbose = option.isVerbose;
    
    if (nargout(objFun) ~= 2) && (nargout(objFun) ~= -1)
        error('objective function should return function value and gradient at x');
    end
    
    try
        [f0, g0] = objFun(x0, true);
    catch err
        error('fail to evaluate the objective function at x0');
    end
    
    if ~(isfinite(f0) && all(isfinite(g0)))
        error('improper initial point x0');
    end
    
    % initialization
    N = numel(x0);
    x = x0; xOld = x0; xEpochOld = x0;
    g = g0; gOld = g0; gEpochOld = g0;
    f = f0; fOld = f0; fEpochOld = f0;
    fDiff = 1e10; xDiff = 1e10;
    numWrongEpoch = 0;
    
    if isempty(option.nonnegMask)
        nonnegMask = false(N, 1);
    else
        nonnegMask = option.nonnegMask;
    end
    
    Mt = zeros(N, 1, clsName); MtEpochOld = Mt;
    Vt = zeros(N, 1, clsName); VtEpochOld = Vt;
    
    m = 1; mOld = 1;
    
    if isVerbose
        fprintf('%d: f = %.3f, df = %.2e, dx = %.2e, eta = %.2e, epoch = %d\n', ...
                    0, fOld, nan, nan, eta, epoch);
    end
    
    while m < maxNumIter
        % calculate momentums
        Mt = b1 .* Mt + (1 - b1) .* gOld;
        Vt = b2 .* Vt + (1 - b2) .* (gOld.^2);
        mth = Mt ./ (1 - b1^m);
        vth = Vt ./ (1 - b2^m);
        
        % learning rate decay based on number of epochs
        if numWrongEpoch >= 1
            decay = 1 / numWrongEpoch^2;
        else
            decay = sqrt(1 - numWrongEpoch);
        end
        eta = option.learningRate * decay;
        
        % calculate directions
        if option.isNesterovAccelerated
            s = eta ./ (sqrt(vth) + epsilon) ...
                .* (b1 .* mth + (1 - b1) .* gOld ./ (1 - b1^m));
        else
            s = eta .* mth ./ (sqrt(vth) + epsilon);
        end
        
        % move to the new point
        x = xOld - s;
        
        % apply nonnegative constraint if any
        mask = nonnegMask & (x < 0);
        x(mask) = 0;
        
%         xDiff = norm(x - xOld, 1) / length(x);
        xDiff = norm(x - xOld) / norm(xOld);
        
        % evaluate f once per epoch, roll back and decay if wrong
        if (m == 1) || (mod(m, epoch) == 0)
            [f, g] = objFun(x, true);
            
            fDiff = abs((f - fOld) / fOld);
            
            if f < fEpochOld                
                % remember evaluation on last epoch
                fEpochOld = f;
                xEpochOld = x;
                gEpochOld = g;
                MtEpochOld = Mt;
                VtEpochOld = Vt;
                mOld = m;
                
                % move x for the correct direction
                xOld = x;
                gOld = g;
                fOld = f;
                
                % decrease #epoch to increase learning rate
                numWrongEpoch = numWrongEpoch - 1;
            else
                if epoch == 1
                    if isVerbose
                        fprintf('wrong step and epoch = 1, do nothing, Mt and Vt will decay\n');
                    end
                    mOld = m;
                else
                    if isVerbose
                        fprintf('wrong step, roll back and reduce epoch\n');
                    end
                    f = fEpochOld;
                    x = xEpochOld;
                    g = gEpochOld;
                    Mt = MtEpochOld;
                    Vt = VtEpochOld;
                    m = mOld;

                    xOld = x;
                    gOld = g;
                    fOld = f;
                    
                    epoch = floor(epoch * 0.8);
                    if epoch < 1
                        epoch = 1;
                    end
                end
                
                % increase #epoch to decrease learning rate
                numWrongEpoch = numWrongEpoch + 1;
            end
        else
            % move x without evaluating f
            [~, g] = objFun(x, false);
            xOld = x;
            gOld = g;
        end
        
        printItv = min([printItv, epoch]);
        if isVerbose
            if (printItv > 0) && ((m <= 1) || (mod(m, printItv) == 0))
                fprintf('%d: f = %.3f, df = %.2e, dx = %.2e, eta = %.2e, epoch = %d\n', ...
                        m, fOld, fDiff, xDiff, eta, epoch);
            end
        end
        
        % check convergence
        if (xDiff < option.tolX) || (fDiff < option.tolFunc)
            if fDiff > option.tolFunc
                [f, ~] = objFun(x, true);
                fDiff = abs((f - fOld) / fOld);
            end
            break;
        end
        
        m = m + 1;
    end
    
    if isVerbose
        if printItv > 0
            if m == maxNumIter
                disp('Reached the max number of iterations');
            else
                disp('Converged');
            end

            disp('Final iteration:');
            fprintf('%d: f = %.3f, df = %.2e, dx = %.2e, eta = %.2e, epoch = %d\n', ...
                        m, fOld, fDiff, xDiff, eta, epoch);
        end
    end
    
    output = struct;
    output.f = f;
    output.numItr = m;
end
