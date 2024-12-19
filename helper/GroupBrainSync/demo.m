clear; close all; clc;

% put the directory of the downloaded folder below or just switch to that folder in Matlab
mainDir = pwd;
addpath(fullfile(mainDir, 'dependency'));

%%
T = 20; V = 50; S = 10;
data = randn(T, V, S);

option = groupBrainSync();
[Y, O] = groupBrainSync(data, option);
