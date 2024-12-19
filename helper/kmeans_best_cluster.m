function [best_parcel, min_sumd_rec] = kmeans_best_cluster(parcels, sumd_rec)
%     to select the best parcellation from multiple k-means trials based on
%     smallest sum of squares error

%     parcels: #vert_num * exp_num
%     sumd_rec: #cluster * exp_num
    
    sumd_rec_sum = sum(sumd_rec);
    [min_sumd_rec, sumd_rec_sum_min_idx] = min(sumd_rec_sum);
    best_parcel = parcels(:, sumd_rec_sum_min_idx);

end