%
% fbn = tsFroNorm(X, mode)
% 
% Description:
%     the Frobenius norm (or square) of the tensor (N-D supported)
% 
% Input:
%     X - the full tensor
%     mode - 1 for the norm and 2 for the squared norm
%
% Output:
%     fbn - the (squared) Frobenius norm
% 
% Copyright:
%     2017-2021 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.0.3
% Date:
%     2021/07/04
%

function fbn = tsFroNorm(X, mode)
    
    if ~exist('mode', 'var') || isempty(mode)
        mode = 1;
    end
    
    X2 = matricize(X, 1);
    
    if mode == 1
        fbn = norm(X2, 'fro');
    elseif mode == 2
        fbn = norm(X2, 'fro') .^ 2;
    end
    
end
