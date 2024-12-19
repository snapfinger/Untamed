function [parcels, sumd_rec] = kmeans_cluster(embed_L, embed_R, cluster_L, cluster_R, km_max_iter, exp_num)
% exp_num: #clustering results from different random initializations

% fsaverage 6 space
    VERT_NUM = 74947;
    parcels = zeros(VERT_NUM, exp_num);
    sumd_rec = zeros(cluster_L + cluster_R, exp_num);
    
    for i = 1: exp_num

        fprintf('%d ', i);

        [km_L, ~, sumd_L] = kmeans(embed_L, cluster_L, 'MaxIter', km_max_iter);
        [km_R, ~, sumd_R] = kmeans(embed_R, cluster_R, 'MaxIter', km_max_iter);
        km_R = km_R + max(km_L);
        
        km = [km_L; km_R];
        parcels(:, i) = km;
        sumd_rec(1: cluster_L, i) = sumd_L;
        sumd_rec(cluster_L+1:end, i) = sumd_R;

        if ~(mod(i, 20))
            fprintf('\n');
        end

    end
    
end