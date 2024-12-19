function embed_rst = netmf_embed_wrapper(graph_adj, negative_sampling, alpha, T, embed_dim)
    
    [un, sn] = netmf_small(graph_adj, T, negative_sampling, embed_dim);
    
    embed_rst = {};
    fprintf('embedding #dim: %d\n', embed_dim);
    embed_cur = un.* (sn .^ alpha);
    embed_rst{1, 1} = embed_dim;
    embed_rst{1, 2} = embed_cur;

end