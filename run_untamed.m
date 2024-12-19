clc; clear;
addpath(genpath('.'));
   
output_folder = './untamed_output/';
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

adj_path = [output_folder 'adj.mat'];
embedding_path = [output_folder 'embed.mat'];
parcel_path = [output_folder 'parcel.mat'];

%% Hyperparameters setup
% see more details in the paper regarding the following 3 hyperparameters

nb=1; % spatial neighborhood constraint
T=7; % NetMF embedding window length
embed_dim=128; % NetMF output embedding dimension

% desired parcel number of left and right hemispheres
cluster_num_left=50;
cluster_num_right=50;

%% graph adjacency matrix construction 

% load NASCAR spatial map
load('artifact/nascar_GSP_rst.mat');

spatial = rst_50{2, 1};
nascar_corr = single(corr(spatial'));

% fsaverage 6 surface space
LEFT_VERT_NUM = 37476;
RIGHT_VERT_NUM = 37471;

% get left and right graph separately
nascar_corr_L = nascar_corr(1: LEFT_VERT_NUM, 1: LEFT_VERT_NUM);
nascar_corr_R = nascar_corr(LEFT_VERT_NUM+1:end, LEFT_VERT_NUM+1:end);

sigma=0.5;
deno = 2 * sigma^2;
nascar_corr_L = exp(nascar_corr_L ./ deno);
nascar_corr_R = exp(nascar_corr_R ./ deno);
nascar_corr_L = (nascar_corr_L + nascar_corr_L') / 2;
nascar_corr_R = (nascar_corr_R + nascar_corr_R') / 2;
clearvars nascar_corr

% add spatial constraint for left and right hemisphere separately
nbinfo = load(['files/fsavg6_nb' num2str(nb) '_index.mat']);
nb_L = nbinfo.nb_L;
nascar_corrL_nb = addSpatialConstraintOneHem(nascar_corr_L, 'L', nb_L);
nb_R = nbinfo.nb_R;
nascar_corrR_nb = addSpatialConstraintOneHem(nascar_corr_R, 'R', nb_R);
clearvars nascar_corr_L nascar_corr_R

save(adj_path, 'nascar_corrL_nb', 'nascar_corrR_nb', '-v7.3');

fprintf('Graph construction done\n');

%% NetMF graph embedding
% to save memory, load and process each hemiphsere adjacency matrix seperately

alpha=0.5;
negative_sampling=1;

fprintf('Start generating embeddings for left hemisphere\n');
embed_L = netmf_embed_wrapper(double(nascar_corrL_nb), negative_sampling, alpha, T, embed_dim);
fprintf('\nDone\n');
clearvars nascar_corrL_nb

fprintf('Start generating embeddings for right hemisphere\n');
embed_R = netmf_embed_wrapper(double(nascar_corrR_nb), negative_sampling, alpha, T, embed_dim);
fprintf('Done\n');
clearvars nascar_corrR_nb
save(embedding_path, 'embed_L','embed_R');

%% Clustering to get final parcellation
km_max_iter=20000;
exp_num=500;

[parcels, sumd_exps] = kmeans_cluster(embed_L{1, 2}, embed_R{1, 2}, ...
                                      cluster_num_left, cluster_num_right, ...
                                      km_max_iter, exp_num);
parcels = int32(parcels);
[parcel, ~] = kmeans_best_cluster(parcels, sumd_exps);
save(parcel_path, 'parcel');

fprintf('\nParcellation generated and saved\n');
