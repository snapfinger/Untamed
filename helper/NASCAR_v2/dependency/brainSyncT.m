%
% [Y2, R] = brainSyncT(X, Y)
% 
% Description:
%     A transposed version of BrainSync for convenient call on data V x T
% 
% Input:
%     X - time series of the reference data (V x T)
%     Y - time series of the subject data (V x T)
% 
% Output:
%     Y2 - syncronized subject data (V x T)
%     R - the orthogonal rotation matrix (V x T)
% 
% Copyright:
%     2018-2021 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.0.1
% Date:
%     2021/10/31
%

function [Y2, R] = brainSyncT(X, Y)

    X = X'; Y = Y';
    [Y2, R] = brainSync(X, Y);
    Y2 = Y2';
    
end
