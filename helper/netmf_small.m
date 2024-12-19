function [u_n, s_n] = netmf_small(Adj, T, neg_sampling, ndim)

    [deepwalk_matrix, ~, ~] = direct_compute_deepwalk_matrix(Adj, T, neg_sampling);
    [u_n, s_n] = svd_deepwalk_matrix(deepwalk_matrix, ndim);

end

%%

function [u, s] = svd_deepwalk_matrix(X, ndim)
% choose to be the max wanted dim
    [u, s, ~] = svds(X, ndim);
    s = diag(s)'; 
end


% direct compute deepwalk matrix
function [Y, L, M] = direct_compute_deepwalk_matrix(A, T, b)
%   params
%     %   A: adjacency matrix
%     %   T: window size
%     %   b: # negative samples
% 
%   return:
%       Y: deepwalk embedding
%       L: symmetric normalized laplacian of A
%       M: deepwalk matrix before taking max & log
    n = size(A, 1);
    vol = sum(A, 'all');
    
    [L, d_rt] = myNormLaplacian(A);

    X = eye(n) - L;
    S = zeros(size(X));
    X_power = eye(n);
    
    for i = 1: T
        fprintf('Compute matrix %d-th power\n', i);
        X_power = X_power * X;
        S = S + X_power;
    end
    
    clearvars X X_power A
    
    S = S * (vol / T / b);
    D_rt_inv = d_rt .^ (-1);
    M = D_rt_inv .* S';
    M = M .* D_rt_inv';
    
    Y = sparse(log(max(M, 1)));
    
end

