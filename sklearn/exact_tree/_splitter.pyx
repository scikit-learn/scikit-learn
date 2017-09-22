#cython: boundscheck=False
#cython: cdivision=True
#cython: warparound=False

from libc.math cimport abs

from ._split_record cimport split_record_copy_to
from ._split_record cimport split_record_reset
from ._split_record cimport split_record_expand_record
from ._split_record cimport split_record_copy_to

from ._stats_node cimport stats_node_copy_to
from ._stats_node cimport stats_node_reset
from ._stats_node cimport stats_node_iadd
from ._stats_node cimport stats_node_isub
from ._stats_node cimport stats_node_copy_to

from ._criterion cimport impurity_improvement

cdef:
    double FEATURE_THRESHOLD = 1e-7
    int TREE_UNDEFINED = -2
    int FEAT_UNKNOWN = -3


cdef void splitter_set_nid(Splitter* splitter, int nid):
    splitter[0].split_record.nid = nid
    splitter[0].original_split_record.nid = nid
    splitter[0].best_split_record.nid = nid


cdef void splitter_init(Splitter* splitter,
                        int feature_idx, int start_idx,
                        SplitRecord* split_record,
                        int min_samples_leaf, double min_weight_leaf):
    splitter[0].feature_idx = feature_idx
    splitter[0].start_idx = start_idx
    splitter[0].prev_idx = start_idx

    split_record_copy_to(split_record, &splitter[0].split_record)
    stats_node_copy_to(&splitter[0].split_record.c_stats,
                       &splitter[0].split_record.r_stats)
    splitter[0].split_record.feature = feature_idx
    splitter[0].split_record.pos = start_idx

    split_record_copy_to(&splitter[0].split_record,
                         &splitter[0].best_split_record)
    split_record_copy_to(&splitter[0].split_record,
                         &splitter[0].original_split_record)

    splitter[0].min_samples_leaf = min_samples_leaf
    splitter[0].min_weight_leaf = min_weight_leaf


cdef void splitter_reset(Splitter* splitter,
                         int feature_idx, int start_idx):
    splitter[0].feature_idx = feature_idx
    splitter[0].start_idx = start_idx
    splitter[0].prev_idx = start_idx

    split_record_copy_to(&splitter[0].original_split_record,
                         &splitter[0].split_record)
    splitter[0].split_record.feature = feature_idx
    splitter[0].split_record.pos = start_idx


cdef void splitter_expand(Splitter* parent_splitter,
                          Splitter* left_splitter,
                          Splitter* right_splitter):
    cdef SplitRecord left_split_record, right_split_record
    split_record_expand_record(&parent_splitter[0].best_split_record,
                               &left_split_record, &right_split_record)
    splitter_init(left_splitter, FEAT_UNKNOWN, TREE_UNDEFINED,
                  &left_split_record,
                  parent_splitter[0].min_samples_leaf,
                  parent_splitter[0].min_weight_leaf)
    splitter_init(right_splitter, FEAT_UNKNOWN, TREE_UNDEFINED,
                  &right_split_record,
                  parent_splitter[0].min_samples_leaf,
                  parent_splitter[0].min_weight_leaf)


cdef void splitter_update_stats(Splitter* splitter, double[::1] y,
                                double[::1] sample_weight,
                                int sample_idx):
    cdef:
        double sum_y = (y[sample_idx] * sample_weight[sample_idx])
        double sum_sq_y = (y[sample_idx] ** 2.0 * sample_weight[sample_idx])
        int n_samples = 1
        double sum_weighted_samples = sample_weight[sample_idx]

    stats_node_reset(&splitter[0].stats_sample, sum_y, sum_sq_y,
                     n_samples, sum_weighted_samples)
    stats_node_iadd(&splitter[0].split_record.l_stats,
                    &splitter[0].stats_sample)
    stats_node_isub(&splitter[0].split_record.r_stats,
                    &splitter[0].stats_sample)


cdef void splitter_node_evaluate_split(Splitter* splitter,
                                       float[:, ::1] X,
                                       double[::1] y,
                                       double[::1] sample_weight,
                                       double sum_total_weighted_samples,
                                       int sample_idx):
    cdef:
        int feat_i = splitter[0].feature_idx
        double diff_samples = abs(X[sample_idx, feat_i] -
                                  X[splitter[0].prev_idx, feat_i])
        bint b_samples_var = diff_samples > FEATURE_THRESHOLD

        int min_samples_leaf = splitter[0].min_samples_leaf
        bint b_n_samples = not(
            splitter[0].split_record.l_stats.n_samples <
            min_samples_leaf or
            splitter[0].split_record.r_stats.n_samples <
            min_samples_leaf)

        double min_weight_leaf = splitter[0].min_weight_leaf
        bint b_weight_samples = not(
            splitter[0].split_record.l_stats.sum_weighted_samples <
            min_weight_leaf or
            splitter[0].split_record.r_stats.sum_weighted_samples <
            min_weight_leaf)
        double current_impurity_improvement = 0.0
    if b_samples_var and b_n_samples and b_weight_samples:
        current_impurity_improvement = impurity_improvement(
            &splitter[0].split_record.c_stats,
            &splitter[0].split_record.l_stats,
            &splitter[0].split_record.r_stats,
            sum_total_weighted_samples)

        if (current_impurity_improvement >
                splitter[0].best_split_record.impurity_improvement):
            split_record_reset(&splitter[0].best_split_record,
                               feature=feat_i,
                               pos=splitter[0].prev_idx,
                               threshold=((X[sample_idx, feat_i] +
                                           X[splitter[0].prev_idx, feat_i]) /
                                          2.0),
                               impurity=splitter[0].split_record.impurity,
                               impurity_improvement=current_impurity_improvement,
                               nid=splitter[0].split_record.nid,
                               c_stats=&splitter[0].split_record.c_stats,
                               l_stats=&splitter[0].split_record.l_stats,
                               r_stats=&splitter[0].split_record.r_stats)

    splitter_update_stats(splitter, y, sample_weight, sample_idx)
    splitter[0].prev_idx = sample_idx


cdef void splitter_copy_to(Splitter* src_splitter,
                           Splitter* dst_splitter):
    dst_splitter[0].feature_idx = src_splitter[0].feature_idx
    dst_splitter[0].start_idx = src_splitter[0].start_idx
    dst_splitter[0].prev_idx = src_splitter[0].prev_idx

    split_record_copy_to(&src_splitter[0].split_record,
                         &dst_splitter[0].split_record)
    split_record_copy_to(&src_splitter[0].best_split_record,
                         &dst_splitter[0].best_split_record)
    split_record_copy_to(&src_splitter[0].original_split_record,
                         &dst_splitter[0].original_split_record)

    stats_node_copy_to(&src_splitter[0].stats_sample,
                       &dst_splitter[0].stats_sample)

    dst_splitter[0].min_samples_leaf = src_splitter[0].min_samples_leaf
    dst_splitter[0].min_weight_leaf = src_splitter[0].min_weight_leaf
