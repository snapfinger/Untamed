%
% parProgressTracker(mode, num)
% 
% Description:
%     Track progress for parallel loop
% 
% Input:
%     mode - 's' - start, 'p' - progress or 'e' - end
%     num - total number of iterations, only needed for mode 's'
% 
% Output:
% 
% Copyright:
%     2016-2021 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian Li (Andrew)
% Revision:
%     2.1.9
% Date:
%     2021/07/03
%

function parProgressTracker(mode, num)

    if strcmpi(mode, 's')
        fprintf(['[' repmat('-', 1, num) ']\n']);
        fprintf('[ \n');
    elseif strcmpi(mode, 'p')
        fprintf('\b\b=>\n');
    elseif strcmpi(mode, 'e')
        fprintf('\b\b]\n');
    end
    
end