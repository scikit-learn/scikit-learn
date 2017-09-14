from ._stats_node cimport StatsNode


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


cdef void split_record_reset(SplitRecord* split_record, int feature,
                             int pos, double threshold,
                             double impurity,
                             double impurity_improvement, int nid,
                             StatsNode* c_stats, StatsNode* l_stats,
                             StatsNode* r_stats)


cdef void split_record_clear(SplitRecord* split_record)


cdef void split_record_expand_record(SplitRecord* split_record,
                                     SplitRecord* left_split_record,
                                     SplitRecord* right_split_record)


cdef void split_record_copy_to(SplitRecord* src_split_record,
                               SplitRecord* dst_split_record)
