from ._stats_node cimport StatsNode


cdef double _impurity_mse(StatsNode stats_node)

cpdef double impurity_mse(double sum_y, double sum_sq_y, int n_samples,
                          double sum_weighted_samples)

cpdef double impurity_improvement(StatsNode c_stats,
                                  StatsNode r_stats,
                                  StatsNode l_stats,
                                  double sum_total_weighted_samples,
                                  criterion)
