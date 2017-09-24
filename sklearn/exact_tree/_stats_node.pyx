cdef void stats_node_copy_to(StatsNode* src_stats_node,
                             StatsNode* dst_stats_node):
    dst_stats_node[0].sum_y = src_stats_node[0].sum_y
    dst_stats_node[0].sum_sq_y = src_stats_node[0].sum_sq_y
    dst_stats_node[0].n_samples = src_stats_node[0].n_samples
    dst_stats_node[0].sum_weighted_samples = src_stats_node[0].sum_weighted_samples
