clc; clear;
code_path = '/hd2/research/Connectivity/code/Parcellation/';
addpath(genpath(code_path));
load('files/setup_params_bobo.mat');

clc

% feat = "rsfc";
feat = "nascar";

for sigma = [0.5]
    for surf_nb = [25, 35]
        
        fprintf('hop %d\n', surf_nb);
        
        % if NASCAR spatial feature as the graph
        if strcmp(feat, "nascar")
            load([data_root 'HCP_3T_Resting/HCP_32k/NASCAR/spatial32k_500sub_35nets.mat']);

            corr_diag = diag(corr(spatial1, spatial2));
            wanted_idx = find(corr_diag > 0.5);
            wanted_idx = setdiff(wanted_idx, 1);
            wanted_idx = setdiff(wanted_idx, 5);
            spatial1 = spatial1(:, wanted_idx);

            % 11: group 1, session 1
            % 21: group 2, session 1
            feat_mat11 = corr(spatial1');
            
        % if RSFC as the graph
        elseif strcmp(feat, 'rsfc')
            load([data_root 'HCP_3T_Resting/HCP_32k/HCP_RSFC_32k_G1S1.mat']);
            feat_mat11 = g1s1_funCon;
            clearvars g1s1_funCon
        else
            error('undefined feature!')
        end
        
        load([code_path 'output/neighbor/32k_nb' num2str(surf_nb) '_index.mat']);

        nb_L = nb_L(~cort32k.idxNaNL, :);
        nb_R = nb_R(~cort32k.idxNaNR, :);
        
        % compute affinity / adjacency matrix
        deno = 2 * sigma^2;
        feat2_mat11 = exp(feat_mat11 ./ deno);
        clearvars feat_mat11

        feat3_mat11 = addSpatialConstraintHem(feat2_mat11, nb_L, nb_R, '32k');
        feat3_mat11 = (feat3_mat11 + feat3_mat11') / 2;
        clearvars feat2_mat11
        
        mat11 = sparse(double(feat3_mat11));
        clearvars feat3_mat11
        
        mat11_L = mat11(1: LEFT_VERT_NUM, 1: LEFT_VERT_NUM);
        mat11_R = mat11(LEFT_VERT_NUM+1: end, LEFT_VERT_NUM+1: end);
        clearvars mat11
        
        if strcmp(feat, "nascar")
            save([code_path 'output/adj_matrix/rsfc_11_hop' num2str(surf_nb)  '_corr_exp_sigma' num2str(sigma) '.mat'], ...
                   'mat11_L', 'mat11_R', '-V7.3');
        elseif strcmp(feat, "rsfc")
            save([code_path 'output/adj_matrix/sub500_n32_v32_np11_hop' num2str(surf_nb)  '_corr_exp_sigma' num2str(sigma) '.mat'], ...
                   'mat11_L', 'mat11_R', '-V7.3');
        end
        clearvars mat11_L mat11_R
               
    end
end