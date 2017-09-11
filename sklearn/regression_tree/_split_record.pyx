#cython: boundscheck=False
#cython: cdivision=True
#cython: warparound=False

from libc.math cimport NAN, INFINITY

from ._stats_node cimport StatsNode
from ._criterion cimport _impurity_mse


cdef class SplitRecord:
    """Container of the parent and children statistics at a node

    Attributes
    ----------
    feature: int,
        The index of the feature to apply the split.

    pos: int,
        The sample index where to apply the split.

    threshold: double,
        The threshold of X[pos, feature] to make the split.

    impurity: double,
        The impurity for that specific split.

    nid: int,
        The id of the current node.

    c_stats: StatsNode,
        Current node statistics.

    l_stats: StatsNode,
        Left node statistics.

    r_stats: StatsNode,
        Right node statistics.
    """

    def __init__(self):
        self.feature = 0
        self.pos = 0
        self.threshold = NAN
        self.impurity = INFINITY
        self.impurity_improvement = -INFINITY
        self.nid = 0

        self.c_stats = StatsNode(0., 0., 0, 0.)
        self.l_stats = StatsNode(0., 0., 0, 0.)
        self.r_stats = StatsNode(0., 0., 0, 0.)

    cpdef void init_stats(self,
                          double c_stats_sum_y, double c_stats_sum_sq_y,
                          int c_stats_n_samples,
                          double c_stats_sum_weighted_samples,
                          double l_stats_sum_y, double l_stats_sum_sq_y,
                          int l_stats_n_samples,
                          double l_stats_sum_weighted_samples,
                          double r_stats_sum_y, double r_stats_sum_sq_y,
                          int r_stats_n_samples,
                          double r_stats_sum_weighted_samples):
        """Init only the statistics of the split record"""

        self.c_stats.sum_y = c_stats_sum_y
        self.c_stats.sum_sq_y = c_stats_sum_sq_y
        self.c_stats.n_samples = c_stats_n_samples
        self.c_stats.sum_weighted_samples = c_stats_sum_weighted_samples

        self.l_stats.sum_y = l_stats_sum_y
        self.l_stats.sum_sq_y = l_stats_sum_sq_y
        self.l_stats.n_samples = l_stats_n_samples
        self.l_stats.sum_weighted_samples = l_stats_sum_weighted_samples

        self.r_stats.sum_y = r_stats_sum_y
        self.r_stats.sum_sq_y = r_stats_sum_sq_y
        self.r_stats.n_samples = r_stats_n_samples
        self.r_stats.sum_weighted_samples = r_stats_sum_weighted_samples

    cpdef void reset(self, int feature, int pos,
                     double threshold, double impurity,
                     double impurity_improvement, int nid,
                     StatsNode c_stats, StatsNode l_stats,
                     StatsNode r_stats):
        """Reset the split record"""
        self.feature = feature
        self.pos = pos
        self.threshold = threshold
        self.impurity = impurity
        self.impurity_improvement = impurity_improvement
        self.nid = nid
        c_stats.copy_to(self.c_stats)
        l_stats.copy_to(self.l_stats)
        r_stats.copy_to(self.r_stats)

    cpdef void clear(self):
        """Clear the split record"""
        self.feature = 0
        self.pos = 0
        self.threshold = NAN
        self.impurity = INFINITY
        self.impurity_improvement = -INFINITY
        self.nid = 0
        self.c_stats.clear()
        self.l_stats.clear()
        self.r_stats.clear()

    cpdef expand_record(self):
        """Create two new records from the left and right stats"""
        # create the left child split record
        left_sr = SplitRecord()
        self.l_stats.copy_to(left_sr.c_stats)
        self.l_stats.copy_to(left_sr.r_stats)
        # FIXME stuck with impurity mse for the moment
        left_sr.impurity = _impurity_mse(left_sr.c_stats)

        # create the right child split record
        right_sr = SplitRecord()
        self.r_stats.copy_to(right_sr.c_stats)
        self.r_stats.copy_to(right_sr.r_stats)
        # FIXME stuck with impurity mse for the moment
        right_sr.impurity = _impurity_mse(right_sr.c_stats)

        return left_sr, right_sr

    # FIXME for debugging purpose
    def __str__(self):
        info = ("feature: {}\n"
                "position: {}\n"
                "threshold: {}\n"
                "impurity: {}\n"
                "impurity improvement: {}\n"
                "node id: {}\n"
                "current stats: {}\n"
                "left stats: {}\n"
                "right stats: {}\n".format(self.feature,
                                           self.pos,
                                           self.threshold,
                                           self.impurity,
                                           self.impurity_improvement,
                                           self.nid,
                                           self.c_stats,
                                           self.l_stats,
                                           self.r_stats))
        return info

    cpdef void copy_to(self, SplitRecord dest_split_record):
        self.c_stats.copy_to(dest_split_record.c_stats)
        self.l_stats.copy_to(dest_split_record.l_stats)
        self.r_stats.copy_to(dest_split_record.r_stats)

        dest_split_record.feature = self.feature
        dest_split_record.pos = self.pos
        dest_split_record.threshold = self.threshold
        dest_split_record.impurity = self.impurity
        dest_split_record.impurity_improvement = self.impurity_improvement
        dest_split_record.nid = self.nid
