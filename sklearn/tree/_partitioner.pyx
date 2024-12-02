"""Partition samples in the construction of a tree.

This module contains the algorithms for moving sample indices to
the left and right child node given a split determined by the
splitting algorithm in `_splitter.pyx`.

Partitioning is done in a way that is efficient for both dense data,
and sparse data stored in a Compressed Sparse Column (CSC) format.
"""
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from cython cimport final
from libc.math cimport isnan, log
from libc.stdlib cimport qsort
from libc.string cimport memcpy

import numpy as np
from scipy.sparse import issparse


# Constant to switch between algorithm non zero value extract algorithm
# in SparsePartitioner
cdef float32_t EXTRACT_NNZ_SWITCH = 0.1

# Allow for 32 bit float comparisons
cdef float32_t INFINITY_32t = np.inf


@final
cdef class DensePartitioner:
    """Partitioner specialized for dense data.

    Note that this partitioner is agnostic to the splitting strategy (best vs. random).
    """
    def __init__(
        self,
        const float32_t[:, :] X,
        intp_t[::1] samples,
        float32_t[::1] feature_values,
        const uint8_t[::1] missing_values_in_feature_mask,
    ):
        self.X = X
        self.samples = samples
        self.feature_values = feature_values
        self.missing_values_in_feature_mask = missing_values_in_feature_mask

    cdef inline void init_node_split(self, intp_t start, intp_t end) noexcept nogil:
        """Initialize splitter at the beginning of node_split."""
        self.start = start
        self.end = end
        self.n_missing = 0

    cdef inline void sort_samples_and_feature_values(
        self, intp_t current_feature
    ) noexcept nogil:
        """Simultaneously sort based on the feature_values.

        Missing values are stored at the end of feature_values.
        The number of missing values observed in feature_values is stored
        in self.n_missing.
        """
        cdef:
            intp_t i, current_end
            float32_t[::1] feature_values = self.feature_values
            const float32_t[:, :] X = self.X
            intp_t[::1] samples = self.samples
            intp_t n_missing = 0
            const uint8_t[::1] missing_values_in_feature_mask = self.missing_values_in_feature_mask

        # Sort samples along that feature; by copying the values into an array and
        # sorting the array in a manner which utilizes the cache more effectively.
        if missing_values_in_feature_mask is not None and missing_values_in_feature_mask[current_feature]:
            i, current_end = self.start, self.end - 1
            # Missing values are placed at the end and do not participate in the sorting.
            while i <= current_end:
                # Finds the right-most value that is not missing so that
                # it can be swapped with missing values at its left.
                if isnan(X[samples[current_end], current_feature]):
                    n_missing += 1
                    current_end -= 1
                    continue

                # X[samples[current_end], current_feature] is a non-missing value
                if isnan(X[samples[i], current_feature]):
                    samples[i], samples[current_end] = samples[current_end], samples[i]
                    n_missing += 1
                    current_end -= 1

                feature_values[i] = X[samples[i], current_feature]
                i += 1
        else:
            # When there are no missing values, we only need to copy the data into
            # feature_values
            for i in range(self.start, self.end):
                feature_values[i] = X[samples[i], current_feature]

        sort(&feature_values[self.start], &samples[self.start], self.end - self.start - n_missing)
        self.n_missing = n_missing

    cdef inline void find_min_max(
        self,
        intp_t current_feature,
        float32_t* min_feature_value_out,
        float32_t* max_feature_value_out,
    ) noexcept nogil:
        """Find the minimum and maximum value for current_feature.

        Missing values are stored at the end of feature_values. The number of missing
        values observed in feature_values is stored in self.n_missing.
        """
        cdef:
            intp_t p, current_end
            float32_t current_feature_value
            const float32_t[:, :] X = self.X
            intp_t[::1] samples = self.samples
            float32_t min_feature_value = INFINITY_32t
            float32_t max_feature_value = -INFINITY_32t
            float32_t[::1] feature_values = self.feature_values
            intp_t n_missing = 0
            const uint8_t[::1] missing_values_in_feature_mask = self.missing_values_in_feature_mask

        # We are copying the values into an array and finding min/max of the array in
        # a manner which utilizes the cache more effectively. We need to also count
        # the number of missing-values there are.
        if missing_values_in_feature_mask is not None and missing_values_in_feature_mask[current_feature]:
            p, current_end = self.start, self.end - 1
            # Missing values are placed at the end and do not participate in the
            # min/max calculation.
            while p <= current_end:
                # Finds the right-most value that is not missing so that
                # it can be swapped with missing values towards its left.
                if isnan(X[samples[current_end], current_feature]):
                    n_missing += 1
                    current_end -= 1
                    continue

                # X[samples[current_end], current_feature] is a non-missing value
                if isnan(X[samples[p], current_feature]):
                    samples[p], samples[current_end] = samples[current_end], samples[p]
                    n_missing += 1
                    current_end -= 1

                current_feature_value = X[samples[p], current_feature]
                feature_values[p] = current_feature_value
                if current_feature_value < min_feature_value:
                    min_feature_value = current_feature_value
                elif current_feature_value > max_feature_value:
                    max_feature_value = current_feature_value
                p += 1
        else:
            min_feature_value = X[samples[self.start], current_feature]
            max_feature_value = min_feature_value

            feature_values[self.start] = min_feature_value
            for p in range(self.start + 1, self.end):
                current_feature_value = X[samples[p], current_feature]
                feature_values[p] = current_feature_value

                if current_feature_value < min_feature_value:
                    min_feature_value = current_feature_value
                elif current_feature_value > max_feature_value:
                    max_feature_value = current_feature_value

        min_feature_value_out[0] = min_feature_value
        max_feature_value_out[0] = max_feature_value
        self.n_missing = n_missing

    cdef inline void next_p(self, intp_t* p_prev, intp_t* p) noexcept nogil:
        """Compute the next p_prev and p for iteratiing over feature values.

        The missing values are not included when iterating through the feature values.
        """
        cdef:
            float32_t[::1] feature_values = self.feature_values
            intp_t end_non_missing = self.end - self.n_missing

        while (
            p[0] + 1 < end_non_missing and
            feature_values[p[0] + 1] <= feature_values[p[0]] + FEATURE_THRESHOLD
        ):
            p[0] += 1

        p_prev[0] = p[0]

        # By adding 1, we have
        # (feature_values[p] >= end) or (feature_values[p] > feature_values[p - 1])
        p[0] += 1

    cdef inline intp_t partition_samples(
        self,
        float64_t current_threshold
    ) noexcept nogil:
        """Partition samples for feature_values at the current_threshold."""
        cdef:
            intp_t p = self.start
            intp_t partition_end = self.end - self.n_missing
            intp_t[::1] samples = self.samples
            float32_t[::1] feature_values = self.feature_values

        while p < partition_end:
            if feature_values[p] <= current_threshold:
                p += 1
            else:
                partition_end -= 1

                feature_values[p], feature_values[partition_end] = (
                    feature_values[partition_end], feature_values[p]
                )
                samples[p], samples[partition_end] = samples[partition_end], samples[p]

        return partition_end

    cdef inline void partition_samples_final(
        self,
        intp_t best_pos,
        float64_t best_threshold,
        intp_t best_feature,
        intp_t best_n_missing,
    ) noexcept nogil:
        """Partition samples for X at the best_threshold and best_feature.

        If missing values are present, this method partitions `samples`
        so that the `best_n_missing` missing values' indices are in the
        right-most end of `samples`, that is `samples[end_non_missing:end]`.
        """
        cdef:
            # Local invariance: start <= p <= partition_end <= end
            intp_t start = self.start
            intp_t p = start
            intp_t end = self.end - 1
            intp_t partition_end = end - best_n_missing
            intp_t[::1] samples = self.samples
            const float32_t[:, :] X = self.X
            float32_t current_value

        if best_n_missing != 0:
            # Move samples with missing values to the end while partitioning the
            # non-missing samples
            while p < partition_end:
                # Keep samples with missing values at the end
                if isnan(X[samples[end], best_feature]):
                    end -= 1
                    continue

                # Swap sample with missing values with the sample at the end
                current_value = X[samples[p], best_feature]
                if isnan(current_value):
                    samples[p], samples[end] = samples[end], samples[p]
                    end -= 1

                    # The swapped sample at the end is always a non-missing value, so
                    # we can continue the algorithm without checking for missingness.
                    current_value = X[samples[p], best_feature]

                # Partition the non-missing samples
                if current_value <= best_threshold:
                    p += 1
                else:
                    samples[p], samples[partition_end] = samples[partition_end], samples[p]
                    partition_end -= 1
        else:
            # Partitioning routine when there are no missing values
            while p < partition_end:
                if X[samples[p], best_feature] <= best_threshold:
                    p += 1
                else:
                    samples[p], samples[partition_end] = samples[partition_end], samples[p]
                    partition_end -= 1


@final
cdef class SparsePartitioner:
    """Partitioner specialized for sparse CSC data.

    Note that this partitioner is agnostic to the splitting strategy (best vs. random).
    """
    def __init__(
        self,
        object X,
        intp_t[::1] samples,
        intp_t n_samples,
        float32_t[::1] feature_values,
        const uint8_t[::1] missing_values_in_feature_mask,
    ):
        if not (issparse(X) and X.format == "csc"):
            raise ValueError("X should be in csc format")

        self.samples = samples
        self.feature_values = feature_values

        # Initialize X
        cdef intp_t n_total_samples = X.shape[0]

        self.X_data = X.data
        self.X_indices = X.indices
        self.X_indptr = X.indptr
        self.n_total_samples = n_total_samples

        # Initialize auxiliary array used to perform split
        self.index_to_samples = np.full(n_total_samples, fill_value=-1, dtype=np.intp)
        self.sorted_samples = np.empty(n_samples, dtype=np.intp)

        cdef intp_t p
        for p in range(n_samples):
            self.index_to_samples[samples[p]] = p

        self.missing_values_in_feature_mask = missing_values_in_feature_mask

    cdef inline void init_node_split(self, intp_t start, intp_t end) noexcept nogil:
        """Initialize splitter at the beginning of node_split."""
        self.start = start
        self.end = end
        self.is_samples_sorted = 0
        self.n_missing = 0

    cdef inline void sort_samples_and_feature_values(
        self,
        intp_t current_feature
    ) noexcept nogil:
        """Simultaneously sort based on the feature_values."""
        cdef:
            float32_t[::1] feature_values = self.feature_values
            intp_t[::1] index_to_samples = self.index_to_samples
            intp_t[::1] samples = self.samples

        self.extract_nnz(current_feature)
        # Sort the positive and negative parts of `feature_values`
        sort(&feature_values[self.start], &samples[self.start], self.end_negative - self.start)
        if self.start_positive < self.end:
            sort(
                &feature_values[self.start_positive],
                &samples[self.start_positive],
                self.end - self.start_positive
            )

        # Update index_to_samples to take into account the sort
        for p in range(self.start, self.end_negative):
            index_to_samples[samples[p]] = p
        for p in range(self.start_positive, self.end):
            index_to_samples[samples[p]] = p

        # Add one or two zeros in feature_values, if there is any
        if self.end_negative < self.start_positive:
            self.start_positive -= 1
            feature_values[self.start_positive] = 0.

            if self.end_negative != self.start_positive:
                feature_values[self.end_negative] = 0.
                self.end_negative += 1

        # XXX: When sparse supports missing values, this should be set to the
        # number of missing values for current_feature
        self.n_missing = 0

    cdef inline void find_min_max(
        self,
        intp_t current_feature,
        float32_t* min_feature_value_out,
        float32_t* max_feature_value_out,
    ) noexcept nogil:
        """Find the minimum and maximum value for current_feature."""
        cdef:
            intp_t p
            float32_t current_feature_value, min_feature_value, max_feature_value
            float32_t[::1] feature_values = self.feature_values

        self.extract_nnz(current_feature)

        if self.end_negative != self.start_positive:
            # There is a zero
            min_feature_value = 0
            max_feature_value = 0
        else:
            min_feature_value = feature_values[self.start]
            max_feature_value = min_feature_value

        # Find min, max in feature_values[start:end_negative]
        for p in range(self.start, self.end_negative):
            current_feature_value = feature_values[p]

            if current_feature_value < min_feature_value:
                min_feature_value = current_feature_value
            elif current_feature_value > max_feature_value:
                max_feature_value = current_feature_value

        # Update min, max given feature_values[start_positive:end]
        for p in range(self.start_positive, self.end):
            current_feature_value = feature_values[p]

            if current_feature_value < min_feature_value:
                min_feature_value = current_feature_value
            elif current_feature_value > max_feature_value:
                max_feature_value = current_feature_value

        min_feature_value_out[0] = min_feature_value
        max_feature_value_out[0] = max_feature_value

    cdef inline void next_p(self, intp_t* p_prev, intp_t* p) noexcept nogil:
        """Compute the next p_prev and p for iteratiing over feature values."""
        cdef:
            intp_t p_next
            float32_t[::1] feature_values = self.feature_values

        if p[0] + 1 != self.end_negative:
            p_next = p[0] + 1
        else:
            p_next = self.start_positive

        while (p_next < self.end and
                feature_values[p_next] <= feature_values[p[0]] + FEATURE_THRESHOLD):
            p[0] = p_next
            if p[0] + 1 != self.end_negative:
                p_next = p[0] + 1
            else:
                p_next = self.start_positive

        p_prev[0] = p[0]
        p[0] = p_next

    cdef inline intp_t partition_samples(
        self,
        float64_t current_threshold
    ) noexcept nogil:
        """Partition samples for feature_values at the current_threshold."""
        return self._partition(current_threshold, self.start_positive)

    cdef inline void partition_samples_final(
        self,
        intp_t best_pos,
        float64_t best_threshold,
        intp_t best_feature,
        intp_t n_missing,
    ) noexcept nogil:
        """Partition samples for X at the best_threshold and best_feature."""
        self.extract_nnz(best_feature)
        self._partition(best_threshold, best_pos)

    cdef inline intp_t _partition(self, float64_t threshold, intp_t zero_pos) noexcept nogil:
        """Partition samples[start:end] based on threshold."""
        cdef:
            intp_t p, partition_end
            intp_t[::1] index_to_samples = self.index_to_samples
            float32_t[::1] feature_values = self.feature_values
            intp_t[::1] samples = self.samples

        if threshold < 0.:
            p = self.start
            partition_end = self.end_negative
        elif threshold > 0.:
            p = self.start_positive
            partition_end = self.end
        else:
            # Data are already split
            return zero_pos

        while p < partition_end:
            if feature_values[p] <= threshold:
                p += 1

            else:
                partition_end -= 1

                feature_values[p], feature_values[partition_end] = (
                    feature_values[partition_end], feature_values[p]
                )
                sparse_swap(index_to_samples, samples, p, partition_end)

        return partition_end

    cdef inline void extract_nnz(self, intp_t feature) noexcept nogil:
        """Extract and partition values for a given feature.

        The extracted values are partitioned between negative values
        feature_values[start:end_negative[0]] and positive values
        feature_values[start_positive[0]:end].
        The samples and index_to_samples are modified according to this
        partition.

        The extraction corresponds to the intersection between the arrays
        X_indices[indptr_start:indptr_end] and samples[start:end].
        This is done efficiently using either an index_to_samples based approach
        or binary search based approach.

        Parameters
        ----------
        feature : intp_t,
            Index of the feature we want to extract non zero value.
        """
        cdef intp_t[::1] samples = self.samples
        cdef float32_t[::1] feature_values = self.feature_values
        cdef intp_t indptr_start = self.X_indptr[feature],
        cdef intp_t indptr_end = self.X_indptr[feature + 1]
        cdef intp_t n_indices = <intp_t>(indptr_end - indptr_start)
        cdef intp_t n_samples = self.end - self.start
        cdef intp_t[::1] index_to_samples = self.index_to_samples
        cdef intp_t[::1] sorted_samples = self.sorted_samples
        cdef const int32_t[::1] X_indices = self.X_indices
        cdef const float32_t[::1] X_data = self.X_data

        # Use binary search if n_samples * log(n_indices) <
        # n_indices and index_to_samples approach otherwise.
        # O(n_samples * log(n_indices)) is the running time of binary
        # search and O(n_indices) is the running time of index_to_samples
        # approach.
        if ((1 - self.is_samples_sorted) * n_samples * log(n_samples) +
                n_samples * log(n_indices) < EXTRACT_NNZ_SWITCH * n_indices):
            extract_nnz_binary_search(X_indices, X_data,
                                      indptr_start, indptr_end,
                                      samples, self.start, self.end,
                                      index_to_samples,
                                      feature_values,
                                      &self.end_negative, &self.start_positive,
                                      sorted_samples, &self.is_samples_sorted)

        # Using an index to samples  technique to extract non zero values
        # index_to_samples is a mapping from X_indices to samples
        else:
            extract_nnz_index_to_samples(X_indices, X_data,
                                         indptr_start, indptr_end,
                                         samples, self.start, self.end,
                                         index_to_samples,
                                         feature_values,
                                         &self.end_negative, &self.start_positive)


cdef int compare_SIZE_t(const void* a, const void* b) noexcept nogil:
    """Comparison function for sort.

    This must return an `int` as it is used by stdlib's qsort, which expects
    an `int` return value.
    """
    return <int>((<intp_t*>a)[0] - (<intp_t*>b)[0])


cdef inline void binary_search(const int32_t[::1] sorted_array,
                               int32_t start, int32_t end,
                               intp_t value, intp_t* index,
                               int32_t* new_start) noexcept nogil:
    """Return the index of value in the sorted array.

    If not found, return -1. new_start is the last pivot + 1
    """
    cdef int32_t pivot
    index[0] = -1
    while start < end:
        pivot = start + (end - start) / 2

        if sorted_array[pivot] == value:
            index[0] = pivot
            start = pivot + 1
            break

        if sorted_array[pivot] < value:
            start = pivot + 1
        else:
            end = pivot
    new_start[0] = start


cdef inline void extract_nnz_index_to_samples(const int32_t[::1] X_indices,
                                              const float32_t[::1] X_data,
                                              int32_t indptr_start,
                                              int32_t indptr_end,
                                              intp_t[::1] samples,
                                              intp_t start,
                                              intp_t end,
                                              intp_t[::1] index_to_samples,
                                              float32_t[::1] feature_values,
                                              intp_t* end_negative,
                                              intp_t* start_positive) noexcept nogil:
    """Extract and partition values for a feature using index_to_samples.

    Complexity is O(indptr_end - indptr_start).
    """
    cdef int32_t k
    cdef intp_t index
    cdef intp_t end_negative_ = start
    cdef intp_t start_positive_ = end

    for k in range(indptr_start, indptr_end):
        if start <= index_to_samples[X_indices[k]] < end:
            if X_data[k] > 0:
                start_positive_ -= 1
                feature_values[start_positive_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, start_positive_)

            elif X_data[k] < 0:
                feature_values[end_negative_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, end_negative_)
                end_negative_ += 1

    # Returned values
    end_negative[0] = end_negative_
    start_positive[0] = start_positive_


cdef inline void extract_nnz_binary_search(const int32_t[::1] X_indices,
                                           const float32_t[::1] X_data,
                                           int32_t indptr_start,
                                           int32_t indptr_end,
                                           intp_t[::1] samples,
                                           intp_t start,
                                           intp_t end,
                                           intp_t[::1] index_to_samples,
                                           float32_t[::1] feature_values,
                                           intp_t* end_negative,
                                           intp_t* start_positive,
                                           intp_t[::1] sorted_samples,
                                           bint* is_samples_sorted) noexcept nogil:
    """Extract and partition values for a given feature using binary search.

    If n_samples = end - start and n_indices = indptr_end - indptr_start,
    the complexity is

        O((1 - is_samples_sorted[0]) * n_samples * log(n_samples) +
          n_samples * log(n_indices)).
    """
    cdef intp_t n_samples

    if not is_samples_sorted[0]:
        n_samples = end - start
        memcpy(&sorted_samples[start], &samples[start],
               n_samples * sizeof(intp_t))
        qsort(&sorted_samples[start], n_samples, sizeof(intp_t),
              compare_SIZE_t)
        is_samples_sorted[0] = 1

    while (indptr_start < indptr_end and
           sorted_samples[start] > X_indices[indptr_start]):
        indptr_start += 1

    while (indptr_start < indptr_end and
           sorted_samples[end - 1] < X_indices[indptr_end - 1]):
        indptr_end -= 1

    cdef intp_t p = start
    cdef intp_t index
    cdef intp_t k
    cdef intp_t end_negative_ = start
    cdef intp_t start_positive_ = end

    while (p < end and indptr_start < indptr_end):
        # Find index of sorted_samples[p] in X_indices
        binary_search(X_indices, indptr_start, indptr_end,
                      sorted_samples[p], &k, &indptr_start)

        if k != -1:
            # If k != -1, we have found a non zero value

            if X_data[k] > 0:
                start_positive_ -= 1
                feature_values[start_positive_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, start_positive_)

            elif X_data[k] < 0:
                feature_values[end_negative_] = X_data[k]
                index = index_to_samples[X_indices[k]]
                sparse_swap(index_to_samples, samples, index, end_negative_)
                end_negative_ += 1
        p += 1

    # Returned values
    end_negative[0] = end_negative_
    start_positive[0] = start_positive_


cdef inline void sparse_swap(intp_t[::1] index_to_samples, intp_t[::1] samples,
                             intp_t pos_1, intp_t pos_2) noexcept nogil:
    """Swap sample pos_1 and pos_2 preserving sparse invariant."""
    samples[pos_1], samples[pos_2] = samples[pos_2], samples[pos_1]
    index_to_samples[samples[pos_1]] = pos_1
    index_to_samples[samples[pos_2]] = pos_2


cdef inline void shift_missing_values_to_left_if_required(
    SplitRecord* best,
    intp_t[::1] samples,
    intp_t end,
) noexcept nogil:
    """Shift missing value sample indices to the left of the split if required.

    Note: this should always be called at the very end because it will
    move samples around, thereby affecting the criterion.
    This affects the computation of the children impurity, which affects
    the computation of the next node.
    """
    cdef intp_t i, p, current_end
    # The partitioner partitions the data such that the missing values are in
    # samples[-n_missing:] for the criterion to consume. If the missing values
    # are going to the right node, then the missing values are already in the
    # correct position. If the missing values go left, then we move the missing
    # values to samples[best.pos:best.pos+n_missing] and update `best.pos`.
    if best.n_missing > 0 and best.missing_go_to_left:
        for p in range(best.n_missing):
            i = best.pos + p
            current_end = end - 1 - p
            samples[i], samples[current_end] = samples[current_end], samples[i]
        best.pos += best.n_missing


# Sort n-element arrays pointed to by feature_values and samples, simultaneously,
# by the values in feature_values. Algorithm: Introsort (Musser, SP&E, 1997).
cdef inline void sort(float32_t* feature_values, intp_t* samples, intp_t n) noexcept nogil:
    if n == 0:
        return
    cdef intp_t maxd = 2 * <intp_t>log(n)
    introsort(feature_values, samples, n, maxd)


cdef inline void swap(float32_t* feature_values, intp_t* samples,
                      intp_t i, intp_t j) noexcept nogil:
    # Helper for sort
    feature_values[i], feature_values[j] = feature_values[j], feature_values[i]
    samples[i], samples[j] = samples[j], samples[i]


cdef inline float32_t median3(float32_t* feature_values, intp_t n) noexcept nogil:
    # Median of three pivot selection, after Bentley and McIlroy (1993).
    # Engineering a sort function. SP&E. Requires 8/3 comparisons on average.
    cdef float32_t a = feature_values[0], b = feature_values[n / 2], c = feature_values[n - 1]
    if a < b:
        if b < c:
            return b
        elif a < c:
            return c
        else:
            return a
    elif b < c:
        if a < c:
            return a
        else:
            return c
    else:
        return b


# Introsort with median of 3 pivot selection and 3-way partition function
# (robust to repeated elements, e.g. lots of zero features).
cdef void introsort(float32_t* feature_values, intp_t *samples,
                    intp_t n, intp_t maxd) noexcept nogil:
    cdef float32_t pivot
    cdef intp_t i, l, r

    while n > 1:
        if maxd <= 0:   # max depth limit exceeded ("gone quadratic")
            heapsort(feature_values, samples, n)
            return
        maxd -= 1

        pivot = median3(feature_values, n)

        # Three-way partition.
        i = l = 0
        r = n
        while i < r:
            if feature_values[i] < pivot:
                swap(feature_values, samples, i, l)
                i += 1
                l += 1
            elif feature_values[i] > pivot:
                r -= 1
                swap(feature_values, samples, i, r)
            else:
                i += 1

        introsort(feature_values, samples, l, maxd)
        feature_values += r
        samples += r
        n -= r


cdef inline void sift_down(float32_t* feature_values, intp_t* samples,
                           intp_t start, intp_t end) noexcept nogil:
    # Restore heap order in feature_values[start:end] by moving the max element to start.
    cdef intp_t child, maxind, root

    root = start
    while True:
        child = root * 2 + 1

        # find max of root, left child, right child
        maxind = root
        if child < end and feature_values[maxind] < feature_values[child]:
            maxind = child
        if child + 1 < end and feature_values[maxind] < feature_values[child + 1]:
            maxind = child + 1

        if maxind == root:
            break
        else:
            swap(feature_values, samples, root, maxind)
            root = maxind


cdef void heapsort(float32_t* feature_values, intp_t* samples, intp_t n) noexcept nogil:
    cdef intp_t start, end

    # heapify
    start = (n - 2) / 2
    end = n
    while True:
        sift_down(feature_values, samples, start, end)
        if start == 0:
            break
        start -= 1

    # sort by shrinking the heap, putting the max element immediately after it
    end = n - 1
    while end > 0:
        swap(feature_values, samples, 0, end)
        sift_down(feature_values, samples, 0, end)
        end = end - 1
