%
% Xn = matricize(X, n)
% 
% Description:
%     matricization along the n-th dimension (N-D supported)
% 
% Input:
%     X - the full tensor
%     n - along which dimension
% 
% Output:
%     Xn - the matricized tensor X (a matrix)
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

function Xn = matricize(X, n)
 
    N = ndims(X);
    if (n <= 0) || (n > N)
        error('mode error');
    end
    
    Y = permute(X, [n 1:n-1 n+1:N]);
    
    sz = size(Y);
    Xn = reshape(Y, sz(1), prod(sz(2:end)));
    
end
