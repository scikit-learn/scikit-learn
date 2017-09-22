from ._split_record cimport SplitRecord

from ._stats_node cimport StatsNode

cdef struct Splitter:
    int feature_idx
    int start_idx
    int prev_idx
    SplitRecord split_record
    SplitRecord original_split_record
    SplitRecord best_split_record
    StatsNode stats_sample
    int min_samples_leaf
    double min_weight_leaf


cdef void splitter_set_nid(Splitter* splitter, int nid)


cdef void splitter_init(Splitter* splitter,
                        int feature_idx, int start_idx,
                        SplitRecord* split_record,
                        int min_samples_leaf, double min_weight_leaf)


cdef void splitter_reset(Splitter* splitter,
                         int feature_idx, int start_idx)


cdef void splitter_expand(Splitter* parent_splitter,
                          Splitter* left_splitter,
                          Splitter* right_splitter)


cdef void splitter_update_stats(Splitter* splitter, double[::1] y,
                                double[::1] sample_weight,
                                int sample_idx)


cdef void splitter_node_evaluate_split(Splitter* splitter,
                                       float[:, ::1] X,
                                       double[::1] y,
                                       double[::1] sample_weight,
                                       double sum_total_weighted_samples,
                                       int sample_idx)


cdef void splitter_copy_to(Splitter* src_splitter,
                           Splitter* dst_splitter)
