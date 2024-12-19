%
% B = nBlockDiag(A, n, isSparse)
% 
% Description:
%     construct block diagonal matrix
% 
% Input:
%     A - the dense block matrix
%     n - how many blocks to repeat
%     isSparse - whether use sparse matrix or not (default is true)
% 
% Output:
%     B - the constructed block diagnoal matrix
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

function B = nBlockDiag(A, n, isSparse)

    if ~exist('isSparse', 'var') || isempty(isSparse)
        isSparse = true;
    end
    
    if isSparse
        B = kron(speye(n, n), A);
    else
        B = kron(eye(n, n), A);
    end
    
end
