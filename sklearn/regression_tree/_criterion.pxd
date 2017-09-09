from ._stats_node cimport StatsNode


cdef float _impurity_mse(StatsNode stats_node)

cpdef float impurity_mse(float sum_y, float sum_sq_y, int n_samples,
                         float sum_weighted_samples)

cpdef float impurity_improvement(StatsNode c_stats,
                                 StatsNode r_stats,
                                 StatsNode l_stats,
                                 float sum_total_weighted_samples,
                                 criterion)
