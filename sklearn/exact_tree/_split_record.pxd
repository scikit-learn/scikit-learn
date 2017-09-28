from libc.math cimport NAN, INFINITY

from ._stats_node cimport StatsNode
from ._stats_node cimport stats_node_copy_to
from ._stats_node cimport stats_node_clear

from ._criterion cimport _impurity_mse


cdef struct SplitRecord:
    int feature
    int pos
    float threshold
    float impurity
    float impurity_improvement
    int nid
    StatsNode c_stats
    StatsNode l_stats
    StatsNode r_stats


cdef void split_record_reset(SplitRecord* split_record, int feature,
                             int pos, float threshold,
                             float impurity,
                             float impurity_improvement, int nid,
                             StatsNode* c_stats, StatsNode* l_stats,
                             StatsNode* r_stats)


cdef inline void split_record_clear(SplitRecord* split_record):
    split_record[0].feature = 0
    split_record[0].pos = 0
    split_record[0].threshold = NAN
    split_record[0].impurity = INFINITY
    split_record[0].impurity_improvement = -INFINITY
    split_record[0].nid = 0

    stats_node_clear(&split_record.c_stats)
    stats_node_clear(&split_record.l_stats)
    stats_node_clear(&split_record.r_stats)


cdef inline void split_record_expand_record(SplitRecord* split_record,
                                            SplitRecord* left_split_record,
                                            SplitRecord* right_split_record):
    split_record_clear(left_split_record)
    stats_node_copy_to(&split_record[0].l_stats, &left_split_record[0].c_stats)
    stats_node_copy_to(&split_record[0].l_stats, &left_split_record[0].r_stats)
    left_split_record[0].impurity = _impurity_mse(
        &left_split_record[0].c_stats)

    split_record_clear(right_split_record)
    stats_node_copy_to(&split_record.r_stats, &right_split_record[0].c_stats)
    stats_node_copy_to(&split_record.r_stats, &right_split_record[0].r_stats)
    right_split_record[0].impurity = _impurity_mse(
        &right_split_record[0].c_stats)


cdef void split_record_copy_to(SplitRecord* src_split_record,
                               SplitRecord* dst_split_record)
