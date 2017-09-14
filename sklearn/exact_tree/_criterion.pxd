#cython: cdivision=True
from ._stats_node cimport StatsNode


cdef inline double _impurity_mse(StatsNode* stats_node):
    cdef double impurity
    impurity = (stats_node[0].sum_sq_y /
                stats_node[0].sum_weighted_samples)
    impurity -= ((stats_node[0].sum_y /
                  stats_node[0].sum_weighted_samples) ** 2.0)

    return impurity


cdef inline double impurity_improvement(StatsNode* c_stats,
                                        StatsNode* l_stats,
                                        StatsNode* r_stats,
                                        double sum_total_weighted_samples):
    # FIXME: only using MSE for the moment
    c_impurity = _impurity_mse(c_stats)
    l_impurity = _impurity_mse(l_stats)
    r_impurity = _impurity_mse(r_stats)

    return ((c_stats[0].sum_weighted_samples /
             sum_total_weighted_samples) *
            (c_impurity -
             (l_stats[0].sum_weighted_samples /
              sum_total_weighted_samples * l_impurity) -
             (r_stats[0].sum_weighted_samples /
              sum_total_weighted_samples * r_impurity)))
