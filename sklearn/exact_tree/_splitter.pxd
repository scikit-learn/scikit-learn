from ._split_record cimport SplitRecord

from ._stats_node cimport StatsNode


cdef struct Splitter:
    float[:, ::1] X
    double[::1] y
    double[::1] sample_weight
    double sum_total_weighted_samples
    int feature_idx
    int start_idx
    int prev_idx
    SplitRecord split_record
    SplitRecord best_split_record
    StatsNode stats_sample
    int min_samples_leaf
    double min_weight_leaf


cdef void splitter_init(Splitter* splitter, float[:, ::1] X, double[::1] y,
                        double[::1] sample_weight,
                        double sum_total_weighted_samples,
                        int feature_idx, int start_idx,
                        SplitRecord* split_record,
                        int min_samples_leaf, double min_weight_leaf)


cdef void splitter_reset(Splitter* splitter, int feature_idx, int start_idx,
                         SplitRecord* split_record)


cdef void splitter_update_stats(Splitter* splitter, int sample_idx)


cdef void splitter_node_evaluate_split(Splitter* splitter, int sample_idx)
