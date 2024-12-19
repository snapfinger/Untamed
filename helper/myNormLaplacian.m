function [L_sym, d_rt] = myNormLaplacian(X)
%   L = I - D^(-1/2) * A * D^(-1/2), where
%   D: degree matrix (diagonal)
%   A: adjacency matrix

    D_diag = sum(X, 1);
    D_diag_inv_sqrt = sqrt(1 ./ D_diag)';
    L_sym = eye(size(X, 1)) - D_diag_inv_sqrt .* X .* D_diag_inv_sqrt';
    d_rt = sqrt(D_diag);
    
end