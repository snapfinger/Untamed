%
% [CCD, G] = CCDWrapper(X, U, lambda)
% 
% Description:
%     CORCONDIA metric wrapper
% 
% Input:
%     X - the full tensor
%     U - the decomposed components
%     lambda - the scales corresponding to the components
% 
% Output:
%     CCD - CORCONDIA metric
% 
% Copyright:
%     2017-2021 (c) LCN & NICC, A. A. Martinos Center, MGH & HMS
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.0.3
% Date:
%     2021/10/31
%

function [CCD, G] = CCDWrapper(X, U, lambda)
    
    lambda = lambda(:);
    U{1} = bsxfun(@times, U{1}, lambda');
    [CCD, G] = corcond(X, U, []);
    
end