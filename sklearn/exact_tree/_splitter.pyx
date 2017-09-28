#cython: boundscheck=False
#cython: cdivision=True
#cython: warparound=False

from ._stats_node cimport stats_node_copy_to


cdef void splitter_init(Splitter* splitter,
                        int feature_idx, int start_idx,
                        SplitRecord* split_record,
                        int min_samples_leaf, float min_weight_leaf):
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


cdef void splitter_copy_to(Splitter* src_splitter,
                           Splitter* dst_splitter):
    dst_splitter[0].feature_idx = src_splitter[0].feature_idx
    dst_splitter[0].start_idx = src_splitter[0].start_idx
    dst_splitter[0].prev_idx = src_splitter[0].prev_idx
    dst_splitter[0].X_prev = src_splitter[0].X_prev

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
