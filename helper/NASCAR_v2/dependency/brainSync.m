%
% [Y2, R] = brainSync(X, Y)
% 
% Description:
%     BrainSync: syncronize the subject fMRI data to the reference fMRI data
% 
% Input:
%     X - time series of the reference data (T x V)
%     Y - time series of the subject data (T x V)
% 
% Output:
%     Y2 - syncronized subject data (T x V)
%     R - the orthogonal rotation matrix (T x T)
% 
% Copyright:
%     2018-2021 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian Li (Andrew), Anand A. Joshi
% Revision:
%     1.0.3
% Date:
%     2021/10/31
%

function [Y2, R] = brainSync(X, Y)

    C = X * Y';
    [U, ~, V] = svd(C);
    R = U * V';
    Y2 = R * Y;
    
end