from libc.math cimport NAN, INFINITY

from ._stats_node cimport StatsNode
from ._stats_node cimport stats_node_reset
from ._stats_node cimport stats_node_clear
from ._stats_node cimport stats_node_copy_to

from ._criterion cimport _impurity_mse


cdef struct SplitRecord:
    int feature
    int pos
    double threshold
    double impurity
    double impurity_improvement
    int nid
    StatsNode c_stats
    StatsNode l_stats
    StatsNode r_stats


# FIXME: no need for this function when the tree builder is in cython as well
cdef inline void split_record_init_stats(SplitRecord split_record,
                                         double c_stats_sum_y,
                                         double c_stats_sum_sq_y,
                                         int c_stats_n_samples,
                                         double c_stats_sum_weighted_samples,
                                         double l_stats_sum_y,
                                         double l_stats_sum_sq_y,
                                         int l_stats_n_samples,
                                         double l_stats_sum_weighted_samples,
                                         double r_stats_sum_y,
                                         double r_stats_sum_sq_y,
                                         int r_stats_n_samples,
                                         double r_stats_sum_weighted_samples):
    stats_node_reset(split_record.c_stats, c_stats_sum_y, c_stats_sum_sq_y,
                     c_stats_n_samples, c_stats_sum_weighted_samples)
    stats_node_reset(split_record.l_stats, l_stats_sum_y, l_stats_sum_sq_y,
                     l_stats_n_samples, l_stats_sum_weighted_samples)
    stats_node_reset(split_record.r_stats, r_stats_sum_y, r_stats_sum_sq_y,
                     r_stats_n_samples, r_stats_sum_weighted_samples)


cdef inline void split_record_reset(SplitRecord split_record, int feature,
                                    int pos, double threshold,
                                    double impurity,
                                    double impurity_improvement, int nid,
                                    StatsNode c_stats, StatsNode l_stats,
                                    StatsNode r_stats):
    split_record.feature = feature
    split_record.pos = pos
    split_record.threshold = threshold
    split_record.impurity = impurity
    split_record.impurity_improvement = impurity_improvement

    stats_node_copy_to(c_stats, split_record.c_stats)
    stats_node_copy_to(l_stats, split_record.l_stats)
    stats_node_copy_to(r_stats, split_record.r_stats)


cdef inline void split_record_clear(SplitRecord split_record):
    split_record.feature = 0
    split_record.pos = 0
    split_record.threshold = NAN
    split_record.impurity = INFINITY
    split_record.impurity_improvement = -INFINITY
    split_record.nid = 0

    stats_node_clear(split_record.c_stats)
    stats_node_clear(split_record.l_stats)
    stats_node_clear(split_record.r_stats)


cpdef inline split_record_expand_record(SplitRecord split_record):
    cdef SplitRecord left_split_record
    cdef SplitRecord right_split_record

    split_record_clear(left_split_record)
    split_record_clear(right_split_record)

    stats_node_copy_to(split_record.l_stats, left_split_record.c_stats)
    stats_node_copy_to(split_record.l_stats, left_split_record.r_stats)
    left_split_record.impurity = _impurity_mse(left_split_record.c_stats)

    stats_node_copy_to(split_record.r_stats, right_split_record.c_stats)
    stats_node_copy_to(split_record.r_stats, right_split_record.r_stats)
    right_split_record.impurity = _impurity_mse(right_split_record.c_stats)

    return left_split_record, right_split_record


cdef inline void split_record_copy_to(SplitRecord src_split_record,
                                      SplitRecord dst_split_record):
    dst_split_record.feature = src_split_record.feature
    dst_split_record.pos = src_split_record.pos
    dst_split_record.threshold = src_split_record.threshold
    dst_split_record.impurity = src_split_record.impurity
    dst_split_record.impurity_improvement = src_split_record.impurity_improvement
    dst_split_record.nid = src_split_record.nid

    stats_node_copy_to(src_split_record.c_stats, dst_split_record.c_stats)
    stats_node_copy_to(src_split_record.l_stats, dst_split_record.l_stats)
    stats_node_copy_to(src_split_record.r_stats, dst_split_record.r_stats)
