function newAdj = addSpatialConstraintHem(adj, nb_L, nb_R, resolution)
%     add spatial constraint to adjacency matrix
%     adj: #observation * #observation
%     nb_L: the array of cells telling nb lists, created from andrew's code
%     nb_R: likewise for right hemisphere
    
    if strcmp(resolution, '11k')
        
        nan_map = load('output/NaN_map.mat');
        left_map = nan_map.left_map;
        right_map = nan_map.right_map;

        VERT_LEFT_NONNAN = 10204;
        VERT_RIGHT_NONNAN = 10208;
        
    elseif strcmp(resolution, '32k')
        
        nan_map = load('output/NaN_map_32k.mat');
        left_map = nan_map.left_map;
        right_map = nan_map.right_map;

        VERT_LEFT_NONNAN = 29696;
        VERT_RIGHT_NONNAN = 29716;
        
    else
        
        error('wrong resolution!');
        
    end
    
    newAdj = zeros(size(adj), class(adj));

    
    for i = 1: VERT_LEFT_NONNAN
        cur_nb_list = left_map(nb_L{i});
        cur_nb_list = cur_nb_list(~isnan(cur_nb_list));
        newAdj(cur_nb_list, cur_nb_list) = adj(cur_nb_list, cur_nb_list);
    end
    
    for i = 1: VERT_RIGHT_NONNAN
        cur_nb_list = right_map(nb_R{i});
        cur_nb_list = cur_nb_list(~isnan(cur_nb_list)) + VERT_LEFT_NONNAN;
        newAdj(cur_nb_list, cur_nb_list) = adj(cur_nb_list, cur_nb_list);
    end
    
end