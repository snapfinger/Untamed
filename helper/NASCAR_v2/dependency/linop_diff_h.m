%
% op = linop_diff_h(sz)
% 
% Description:
%     construct a linear differential operator
% 
% Input:
%     sz - size of the matrix
% 
% Output:
%     op - the constructed operator
% 
% Copyright:
%     2017-2021 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.0.1
% Date:
%     2021/07/04
%

function op = linop_diff_h(sz)
    
    n1 = sz(1);
    n2 = sz(2);

    mat = @(x) reshape(x, n1, n2);
    D = @(X) [diff(X,1,2), zeros(n1,1)];
    Dt = @(X) [zeros(n1,1), X(:,1:end-1)] - [X(:,1:end-1), zeros(n1,1)];
    
    TV  = @(x) vec(D(mat(x)));
    TVt = @(x) vec(Dt(mat(x)));

    szW = n1 * n2;
    szW = [szW, szW];
    op = @(x, mode) linop_diff_h_impl(szW, TV, TVt, x, mode);
    
end
    
function y = linop_diff_h_impl(sz, TV, TVt, x, mode)
    switch mode
        case 0, y = sz;
        case 1, y = TV(x);
        case 2, y = TVt(x);
    end
end