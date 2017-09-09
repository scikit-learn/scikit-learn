from ._stats_node cimport StatsNode
from ._criterion cimport impurity_mse


cdef class SplitRecord:

    cdef:
        public int feature
        public int pos
        public double threshold
        public double impurity
        public double impurity_improvement
        public int nid
        public StatsNode c_stats
        public StatsNode l_stats
        public StatsNode r_stats

    cpdef void init_stats(self,
                          double c_stats_sum_y, double c_stats_sum_sq_y,
                          int c_stats_n_samples,
                          double c_stats_sum_weighted_samples,
                          double l_stats_sum_y, double l_stats_sum_sq_y,
                          int l_stats_n_samples,
                          double l_stats_sum_weighted_samples,
                          double r_stats_sum_y, double r_stats_sum_sq_y,
                          int r_stats_n_samples,
                          double r_stats_sum_weighted_samples)

    cpdef void reset(self, int feature, int pos,
                     double threshold, double impurity,
                     double impurity_improvement, int nid,
                     StatsNode c_stats, StatsNode l_stats,
                     StatsNode r_stats)

    cpdef void clear(self)

    cpdef expand_record(self)

    cpdef void copy_to(self, SplitRecord dest_split_record)
