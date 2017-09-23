#cython: boundscheck=False
#cython: cdivision=True
#cython: warparound=False

# from libc.math cimport abs

from ._split_record cimport split_record_expand_record

from ._stats_node cimport stats_node_copy_to
from ._stats_node cimport stats_node_copy_to


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
