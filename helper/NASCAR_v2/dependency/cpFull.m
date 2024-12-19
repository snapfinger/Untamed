%
% Y = cpFull(U, lambda, limitedMem)
% 
% Description:
%     reconstruct the full tensor based on the components (N-D supported)
% 
% Input:
%     U - the decomposed components
%     lambda - the scales corresponding to the components
%     limitedMem - whether use limited memory mode or not 
%
% Output:
%     Y - the reconstructed full tensor
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

function Y = cpFull(U, lambda, limitedMem)
    
    if ~exist('limitedMem', 'var') || isempty(limitedMem)
        limitedMem = true;
    end
    
    N = length(U);
    sz = zeros(1, N);
    R = size(U{1}, 2);
    clsName = class(U{1});
    
    if ~exist('lambda', 'var') || isempty(lambda)
        lambda = ones(R, 1, clsName);
    end
    
    for m = 1:N
        sz(m) = size(U{m}, 1);
    end
    
    if limitedMem
        Y = zeros(sz, clsName);
        for m = 1:R
            T = U{1}(:, m);
            for n = 2:N
                s = ones(1, N);
                s(n) = sz(n);
                v = reshape(U{n}(:, m), s);
                T = T .* v;
            end
            Y = Y + lambda(m) .* T;
        end
        
    else
        Y = reshape(krProd(U(end:-1:1)) * lambda, sz);
    end
end
