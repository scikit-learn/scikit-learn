cdef class StatsNode:

    cdef:
        readonly double sum_y, sum_sq_y, sum_weighted_samples
        readonly int n_samples

    cdef void reset(self, double sum_y, double sum_sq_y, int n_samples,
                     double sum_weighted_samples)

    cdef void clear(self)

    cdef void copy_to(self, StatsNode dest_stats_node)
