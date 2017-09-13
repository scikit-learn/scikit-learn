cdef struct StatsNode:
    double sum_y
    double sum_sq_y
    int n_samples
    double sum_weighted_samples


cdef StatsNode stats_node_reset(StatsNode stats_node,
                                double sum_y, double sum_sq_y, int n_samples,
                                double sum_weighted_samples)


cdef StatsNode stats_node_clear(StatsNode stats_node)


cdef StatsNode stats_node_iadd(StatsNode l_stats_node,
                          StatsNode r_stats_node)


cdef StatsNode stats_node_isub(StatsNode l_stats_node,
                          StatsNode r_stats_node)


cdef StatsNode stats_node_copy_to(StatsNode src_stats_node,
                                  StatsNode dst_stats_node)
