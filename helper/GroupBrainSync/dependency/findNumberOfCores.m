%
% numCores = findNumberOfCores()
% 
% Description:
%     find number of available cpu cores
% 
% Input:
% 
% Output:
%     numCores - number of cpu cores
% 
% Copyright:
%     2014-2021 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian Li (Andrew)
% Revision:
%     1.0.3
% Date:
%     2018/01/22
%

function numCores = findNumberOfCores()
    numCores = feature('numCores');
end