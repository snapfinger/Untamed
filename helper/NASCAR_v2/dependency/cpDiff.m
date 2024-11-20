%
% [df, relerr, varExp] = cpDiff(X, U, lambda)
% 
% Description:
%     calculate the difference between tensor and approximations (N-D supported)
% 
% Input:
%     X - the full tensor
%     U - the decomposed components
%     lambda - the scales corresponding to the components
% 
% Output:
%     df - the Frobenius norm difference
%     relerr - relative error to the norm of the original tensor X
%     varExp - variance explained
% 
% Copyright:
%     2017-2021 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.0.7
% Date:
%     2021/07/04
%

function [df, relerr, varExp] = cpDiff(X, U, lambda)

    Xnorm = tsFroNorm(X);
    
    A = matricize(X, 1) - U{1} * diag(lambda) * krProd(U(end:-1:2))';
    df = norm(A, 'fro');
        
    relerr = df ./ Xnorm;
    varExp = 1 - relerr.^2;
    
end
