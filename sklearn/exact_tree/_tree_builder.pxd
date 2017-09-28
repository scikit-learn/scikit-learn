from ._tree cimport Tree
from ._tree cimport TreeBuilder

from ._splitter cimport Splitter


cdef void weighted_sum_y(float[::1] y, float[::1] sample_weight,
                         float* p_sum_y, float* p_sum_sq_y)


cdef class ExactTreeBuilder(TreeBuilder):

    cdef:
        int min_samples_split
        int min_samples_leaf
        float min_weight_leaf
        float min_impurity_split
        int max_depth
        int max_features

    cpdef build(self, Tree tree, float[:, ::1] X, int[::1, :] X_idx_sorted,
                float[::1] y, float[::1] sample_weight,
                float sum_total_weighted_samples)

    cdef inline bint _check_minimum_stats(self, Splitter* splitter):
        cdef:
            b_impurity = (
                splitter[0].split_record.impurity >
                self.min_impurity_split)
            b_samples_split = (
                splitter[0].split_record.c_stats.n_samples >=
                self.min_samples_split)
            b_samples_leaf = (
                splitter[0].split_record.c_stats.n_samples >=
                self.min_samples_leaf)

        return b_impurity and b_samples_leaf and b_samples_split
