cdef void stats_node_clear(StatsNode* stats_node):
    stats_node[0].sum_y = 0.0
    stats_node[0].sum_sq_y = 0.0
    stats_node[0].n_samples = 0
    stats_node[0].sum_weighted_samples = 0.0


cdef void stats_node_copy_to(StatsNode* src_stats_node,
                             StatsNode* dst_stats_node):
    dst_stats_node[0].sum_y = src_stats_node[0].sum_y
    dst_stats_node[0].sum_sq_y = src_stats_node[0].sum_sq_y
    dst_stats_node[0].n_samples = src_stats_node[0].n_samples
    dst_stats_node[0].sum_weighted_samples = src_stats_node[0].sum_weighted_samples
