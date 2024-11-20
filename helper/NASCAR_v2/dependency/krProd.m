%
% KR = krProd(U)
% 
% Description:
%     implementation of the Khatri-Rao product (N-D supported)
% 
% Input:
%     U - the components in cell structure
% 
% Output:
%     KR - the product
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

function KR = krProd(U)
    
    N = length(U);
    cols = zeros(N, 1);
    for m = 1:N
        cols(m) = size(U{m}, 2);
    end
    
    if(~all(cols == cols(1)))
        error('number of columns does not match');
    end
    
    KR = U{1};
    R = cols(1);
    
    for m = 2:N
        KR = bsxfun(@times, reshape(U{m}, [], 1, R), reshape(KR, 1, [], R));
    end
    
    KR = reshape(KR, [], R);
    
end
