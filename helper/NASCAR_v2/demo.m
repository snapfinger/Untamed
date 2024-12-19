clear; close all; clc;

% put the directory of the downloaded folder below or just switch to that folder in Matlab
mainDir = pwd;
addpath(fullfile(mainDir, 'dependency'));

% (optional) uncomment the following two lines and provide the TFOCS directory
% tfocsDir = '';
% addpath(tfocsDir);

% (optional) uncomment the following two lines and provide the nway toolbox directory
% nwayDir = '';
% addpath(nwayDir);

%% simulate data
N = 1000; T = 500; S = 20;
R = 5;
UGT = {randn(N, R); randn(T, R); rand(S, R)};
lambdaGT = [R:-1:1] * 10;
data = cpFull(UGT, lambdaGT);

ns = randn(size(data)) * 1;
data = data + ns;

%% BrainSync temporal alignment
for m = 2:S
    data(:, :, m) = brainSyncT(data(:, :, 1), data(:, :, m));
    m
end

%% NASCAR algorithm (use opt inside srscpd framework)
option = srscpd('opt');
option.nonnegative = [0 0 1];
option.rankOneOptALS.useTFOCS = false; % if use TFOCS or not to solve rank 1 problem
option.optAlg.normRegParam = 0.001; % regularization parameter on the norm of each mode
option.optAlg.optSolver.printItv = 10; % print interval
option.optAlg.optSolver.learningRate = 1e-3; % Nadam learning rate
option.optAlg.optSolver.maxNumIter = 2000; % max number of iterations in Nadam
option.optAlg.cacheMTS = true; % turn this on if you have enough memory to speed up
% save result at the end of each iteration during the decomposition
% option.saveToFile = 'file_to_save_result';

% save logs (what printed in the command window) to a file
% option.logFile = 'log_file';

% if NASCAR paused or is interrupted for some reason, it can be resumed
% by providing the results from lower ranks
% option.resumeFrom = 'lower_rank_results';

result = srscpd(data, R, option);

%% the decomposition results are stored in U
U2 = result(R).U;
lambda2 = result(R).Lambda;

c = corr(UGT{1}, U2{1});
figure, imagesc(abs(c));

% the sign of the spatial and temporal mode is arbitrary, flip them to
% match for simulation
for m = 1:R
    if c(m, m) < 0
        U2{1}(:, m) = -U2{1}(:, m);
        U2{2}(:, m) = -U2{2}(:, m);
    end
end

%% plot to check the results
mode = 1;
r = 1;
figure;
plot(UGT{mode}(:, r) * lambdaGT(r).^(1/3));
hold on; grid on;
plot(U2{mode}(:, r) * lambda2(r).^(1/3));

disp(corr(UGT{mode}(:, r), U2{mode}(:, r)));
