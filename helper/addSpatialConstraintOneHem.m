function newAdj = addSpatialConstraintOneHem(adj, hem, nb)
%     add spatial constraint to adjacency matrix
%     adj: #observation * #observation
%     hem: 'L' or 'R'
%     nb: the array of cells telling nb lists of intended hemispehre, 
%         created from andrew's code
    
    nan_map = load('files/NaN_map_fsavg6.mat');
    if strcmp(hem, 'L')
        map = nan_map.mapL_fsavg6;
        VERT_NONNAN = 37476;
    elseif strcmp(hem, 'R')
        map = nan_map.mapR_fsavg6;
        VERT_NONNAN = 37471;
    else
        error('wrong hemisphere input');
    end

    newAdj = zeros(size(adj), class(adj));

    for i = 1: VERT_NONNAN
        cur_nb_list = map(nb{i});
        cur_nb_list = cur_nb_list(~isnan(cur_nb_list));
        newAdj(cur_nb_list, cur_nb_list) = adj(cur_nb_list, cur_nb_list);
    end
    
end