clc; clear; addpath(genpath('.'));
addpath /hd2/sw1/MATLAB18/gifti
addpath('/hd2/sw1/cifti-matlab');
addpath(genpath('/hd2/sw1/bfp/'))

%%
LEFT_VERT_NUM = 29696;
VERT_NUM = 59412;

data_root = '/hd2/research/Connectivity/data/';
code_path = '/hd2/research/Connectivity/code/Parcellation/';
rst_root = [code_path 'output/adj_matrix/'];

%% embedding
exp_num = 50;
negative_sampling = 1; sigma = 0.5;
                
for sess = [11]
    fprintf('sess: %d\n', sess);
    
    for alpha = [0.5]
        for surf_nb = [7]
            
            fprintf('hop: %d\n', surf_nb);
            
            for T = [7]
                
                feat_path = [rst_root 'sub500_n32_v32_np11_hop' num2str(surf_nb) ...
                             '_corr_exp_sigma' num2str(sigma) '.mat'];
                load(feat_path);
                
                for embed_dim = [200]
                    
                    fprintf('embedding #dim: %d\n', embed_dim);

                    if sess == 11
                        dp_eb_L = netmf_small(mat11_L, T, negative_sampling, embed_dim, alpha);
                        dp_eb_R = netmf_small(mat11_R, T, negative_sampling, embed_dim, alpha);
                    elseif sess == 21
                        dp_eb_L = netmf_small(mat21_L, T, negative_sampling, embed_dim, alpha);
                        dp_eb_R = netmf_small(mat21_R, T, negative_sampling, embed_dim, alpha);
                    end
                    
                    save([rst_root 'embed/netmf_' num2str(sess) '_hop' num2str(surf_nb) '_sig' num2str(sigma) ...
                          '_d' num2str(embed_dim)  ...
                          '_T' num2str(T) '_alpha' num2str(alpha) '.mat'], 'dp_eb_L', 'dp_eb_R');
                end

            end
        end
    end
end