from ._tree cimport Tree
from ._tree cimport TreeBuilder

from ._splitter cimport Splitter


cdef void weighted_sum_y(double[::1] y, double[::1] sample_weight,
                         double* p_sum_y, double* p_sum_sq_y)


cdef class ExactTreeBuilder(TreeBuilder):

    cdef:
        int min_samples_split
        int min_samples_leaf
        double min_weight_leaf
        double min_impurity_split
        int max_depth
        int max_features

    cpdef build(self, Tree tree, float[:, ::1] X, int[::1, :] X_idx_sorted,
                double[::1] y, double[::1] sample_weight,
                double sum_total_weighted_samples)

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
