%
% [err, relerr] = cpDiffX(U1, U2, base)
% 
% Description:
%     calculate norm difference between solutions of tensor decomposition
% 
% Input:
%     U1 - tensor solution 1
%     U2 - tensor solution 2
%     base - which one used as the base to calculate the relative norm
%
% Output:
%     err - the absolute norm differnece
%     relerr - relative norm difference
% 
% Copyright:
%     2019-2021 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.0.2
% Date:
%     2021/07/04
%

function [err, relerr] = cpDiffX(U1, U2, base)

    if ~exist('base', 'var') || isempty(base)
        base = 2;
    end        
    
    N = length(U1);
    clsName = class(U1{1});
    
    xDiffs = zeros(N, 1, clsName);
    oldNorms = zeros(N, 1, clsName);
    for n = 1:N
        xDiffs(n) = norm(U1{n} - U2{n}, 'fro');
        if base == 1
            oldNorms(n) = norm(U1{n}, 'fro');
        elseif base == 2
            oldNorms(n) = norm(U2{n}, 'fro');
        end
    end
    
    err = norm(xDiffs);
    relerr = err / norm(oldNorms);

end