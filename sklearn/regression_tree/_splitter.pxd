from sklearn.regression_tree._stats_node cimport StatsNode
from sklearn.regression_tree._split_record cimport SplitRecord

cdef class Splitter:

    cdef:
        float[:, ::1] X
        double[::1] y, sample_weight
        double sum_total_weighted_samples

        int feature_idx, start_idx, prev_idx

        int min_samples_leaf
        double min_weight_leaf

        SplitRecord split_record
        readonly SplitRecord best_split_record

        StatsNode stats_samples

    cpdef void reset(self, int feature_idx, int start_idx,
                     SplitRecord split_record)

    cpdef void update_stats(self, int sample_idx)

    cpdef void node_evaluate_split(self, int sample_idx)
