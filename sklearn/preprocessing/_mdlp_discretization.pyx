# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3


# Author: Juan Carlos Alfaro Jiménez <JuanCarlos.Alfaro@uclm.es>
# License: BSD


import numpy as np
from libc.math cimport log, pow
from libc.math cimport INFINITY
from libc.stdlib cimport calloc, free, realloc
from libc.string cimport memcpy, memset
from numpy import float64 as DTYPE
from numpy import int32 as ITYPE
from numpy import intp as SIZE


# =============================================================================
# Constants
# =============================================================================

# Initial size for the stack
cdef SIZE_t INITIAL_STACK_SIZE = 10

# Maximum bin size
cdef SIZE_t MAX_BINS = 255


# =============================================================================
# Functions
# =============================================================================

def find_binning_thresholds(const DTYPE_t[:, :] X,
                            const SIZE_t[:, :] X_idx_sorted,
                            const ITYPE_t[:, :] y,
                            SIZE_t n_outputs,
                            SIZE_t[:] n_classes):
    """Find the thresholds used to discretize the data."""
    cdef SIZE_t n_samples = X.shape[0]
    cdef SIZE_t n_features = X.shape[1]
    cdef SIZE_t n_elements

    cdef DTYPE_t* sum_total = NULL
    cdef SIZE_t sum_stride = 0

    cdef np.ndarray[SIZE_t, ndim=1] n_bins
    cdef np.ndarray[DTYPE_t, ndim=2] bin_edges
    cdef np.ndarray[object, ndim=1] bin_edges_output

    cdef SIZE_t label
    cdef SIZE_t output
    cdef SIZE_t feature

    n_bins = np.ones(n_features, dtype=SIZE)

    # Need to use a maximum number of bins to initialize the view
    bin_edges = np.zeros((n_features, MAX_BINS), dtype=DTYPE)
    bin_edges_output = np.zeros(n_features, dtype=object)

    # Compute the maximal stride of all targets
    for output in range(n_outputs):
        if n_classes[output] > sum_stride:
            sum_stride = n_classes[output]

    n_elements = n_outputs * sum_stride
    safe_calloc(&sum_total, n_elements)

    # Compute the sum of the count of each label to cache
    # the results and avoid intermediate computations
    for sample in range(n_samples):
        for output in range(n_outputs):
            label = y[sample, output]
            sum_total[output * sum_stride + label] += 1.0

    for feature in range(n_features):
        _find_binning_thresholds(X, feature, X_idx_sorted,
                                 y, sum_total, sum_stride,
                                 n_outputs, n_classes,
                                 bin_edges, n_bins)

        # The cut points must be sorted because of the stack
        bin_edges_output[feature] = np.sort(
            bin_edges[feature, :n_bins[feature] - 1])

    free(sum_total)

    return bin_edges_output, n_bins


cdef SIZE_t _find_binning_thresholds(const DTYPE_t[:, :] X,
                                     SIZE_t feature,
                                     const SIZE_t[:, :] X_idx_sorted,
                                     const ITYPE_t[:, :] y,
                                     DTYPE_t* sum_total,
                                     SIZE_t sum_stride,
                                     SIZE_t n_outputs,
                                     const SIZE_t[:] n_classes,
                                     DTYPE_t[:, :] bin_edges,
                                     SIZE_t[:] n_bins):
    """Find the thresholds used to discretize the feature data."""
    cdef SIZE_t n_samples = X.shape[0]

    cdef Stack stack = Stack(INITIAL_STACK_SIZE)
    cdef StackRecord stack_record

    cdef SIZE_t start
    cdef SIZE_t end
    cdef SIZE_t cut_point
    cdef SIZE_t index = 0

    cdef DTYPE_t* sum_left = NULL
    cdef DTYPE_t* sum_right = NULL
    cdef SIZE_t n_elements = n_outputs * sum_stride

    with nogil:
        # Push the root record onto the stack
        stack.push(0, n_samples, sum_total, sum_stride, n_outputs)

        safe_calloc(&sum_left, n_elements)
        safe_calloc(&sum_right, n_elements)
        safe_calloc(&stack_record.sum_total, n_elements)

        # Recursive partition
        while stack.top > 0:
            # Pop the top record from the stack
            stack.pop(&stack_record, sum_stride, n_outputs)

            start = stack_record.start
            end = stack_record.end
            sum_total = stack_record.sum_total

            cut_point = find_cut_point(X, feature, X_idx_sorted,
                                       y, sum_total, &sum_left,
                                       &sum_right, sum_stride,
                                       n_outputs, n_classes,
                                       start, end)
            if not cut_point:
                # No good cut point
                continue

            accept_split = is_split_accepted(sum_total, sum_left,
                                             sum_right, sum_stride,
                                             n_outputs, n_classes,
                                             start, cut_point, end)

            if not accept_split:
                # No split accepted
                continue

            # Push the left and right part onto the stack
            stack.push(start, cut_point, sum_left, sum_stride, n_outputs)
            stack.push(cut_point, end, sum_right, sum_stride, n_outputs)

            # Set the threshold for the cut point and increase the pointer
            bin_edges[feature, n_bins[feature] - 1] = get_threshold(
                X, feature, X_idx_sorted, cut_point)

            n_bins[feature] += 1

        free(stack_record.sum_total)
        free(sum_left)
        free(sum_right)


cdef SIZE_t find_cut_point(const DTYPE_t[:, :] X,
                           SIZE_t feature,
                           const SIZE_t[:, :] X_idx_sorted,
                           const ITYPE_t[:, :] y,
                           DTYPE_t* sum_total,
                           DTYPE_t** sum_left,
                           DTYPE_t** sum_right,
                           SIZE_t sum_stride,
                           SIZE_t n_outputs,
                           const SIZE_t[:] n_classes,
                           SIZE_t start,
                           SIZE_t end) nogil except -1:
    """Find the cut point between the start and end points."""
    cdef DTYPE_t n_samples = end - start
    cdef DTYPE_t n_samples_left
    cdef DTYPE_t n_samples_right

    cdef SIZE_t sample
    cdef SIZE_t output
    cdef SIZE_t label
    cdef SIZE_t pos = start

    cdef SIZE_t label_index
    cdef SIZE_t prev_sample_index
    cdef SIZE_t curr_sample_index
    cdef DTYPE_t prev_sample_value
    cdef DTYPE_t curr_sample_value

    cdef DTYPE_t entropy
    cdef DTYPE_t entropy_left
    cdef DTYPE_t entropy_right
    cdef DTYPE_t best_entropy = INFINITY

    cdef DTYPE_t* best_sum_left = NULL
    cdef DTYPE_t* best_sum_right = NULL

    cdef SIZE_t n_elements = n_outputs * sum_stride
    cdef SIZE_t n_bytes = n_elements * sizeof(DTYPE_t)

    cdef DTYPE_t count_label
    cdef SIZE_t cut_point = 0  # 0 for no cut point

    # Reset the left and right part statistics
    memset(sum_left[0], 0, n_bytes)
    memset(sum_right[0], 0, n_bytes)

    # Allocate memory for the best left and right part
    # statistics to quite optimize the binning process
    safe_calloc(&best_sum_left, n_elements)
    safe_calloc(&best_sum_right, n_elements)

    for sample in range(start + 1, end):
        prev_sample_index = X_idx_sorted[sample - 1, feature]
        curr_sample_index = X_idx_sorted[sample, feature]
        prev_sample_value = X[prev_sample_index, feature]
        curr_sample_value = X[curr_sample_index, feature]

        # Do not cut samples with the same values
        if prev_sample_value == curr_sample_value:
            continue

        update_statistics(feature, X_idx_sorted, y,
                          sum_total, sum_left, sum_right,
                          sum_stride, n_outputs, n_classes,
                          pos, sample)

        # Update the position
        pos = sample

        # Reset the left and right part statistics
        entropy_left = 0.0
        entropy_right = 0.0

        # Update the number of samples in the left and right part
        n_samples_left = pos - start
        n_samples_right = n_samples - n_samples_left

        for output in range(n_outputs):
            for label in range(n_classes[output]):
                label_index = output * sum_stride + label

                # Compute the entropy for the left part
                count_label = sum_left[0][label_index] / n_samples_left

                if count_label > 0.0:
                    entropy_left -= count_label * log2(count_label)

                # Compute the entropy for the right part
                count_label = sum_right[0][label_index] / n_samples_right

                if count_label > 0.0:
                    entropy_right -= count_label * log2(count_label)

        entropy_left /= n_outputs
        entropy_right /= n_outputs

        entropy = (n_samples_left / n_samples * entropy_left
                   + n_samples_right / n_samples * entropy_right)

        if entropy <= best_entropy:
            best_entropy = entropy
            cut_point = sample

            # Store the best left and right part statistics
            memcpy(best_sum_left, sum_left[0], n_bytes)
            memcpy(best_sum_right, sum_right[0], n_bytes)

    # Copy the best left and right part statistics
    memcpy(sum_left[0], best_sum_left, n_bytes)
    memcpy(sum_right[0], best_sum_right, n_bytes)

    return cut_point


cdef void update_statistics(SIZE_t feature,
                            const SIZE_t[:, :] X_idx_sorted,
                            const ITYPE_t[:, :] y,
                            DTYPE_t* sum_total,
                            DTYPE_t** sum_left,
                            DTYPE_t** sum_right,
                            SIZE_t sum_stride,
                            SIZE_t n_outputs,
                            const SIZE_t[:] n_classes,
                            SIZE_t pos,
                            SIZE_t new_pos) nogil:
    """Update the left and right part statistics."""
    cdef SIZE_t sample
    cdef SIZE_t output
    cdef SIZE_t label

    cdef SIZE_t sample_index
    cdef SIZE_t label_index

    for sample in range(pos, new_pos):
        sample_index = X_idx_sorted[sample, feature]

        # Update the left part statistics
        for output in range(n_outputs):
            label_index = (output * sum_stride
                           + y[sample_index, output])

            sum_left[0][label_index] += 1

    # Update the right part statistics
    for output in range(n_outputs):
        for label in range(n_classes[output]):
            label_index = output * sum_stride + label

            sum_right[0][label] = (sum_total[label_index]
                                   - sum_left[0][label_index])


cdef BOOL_t is_split_accepted(DTYPE_t* sum_total,
                              DTYPE_t* sum_left,
                              DTYPE_t* sum_right,
                              SIZE_t sum_stride,
                              SIZE_t n_outputs,
                              const SIZE_t[:] n_classes,
                              SIZE_t start,
                              SIZE_t pos,
                              SIZE_t end) nogil:
    """Check if the split is accepted."""
    cdef DTYPE_t n_samples_node = end - start
    cdef DTYPE_t n_samples_left = pos - start
    cdef DTYPE_t n_samples_right = n_samples_node - n_samples_left

    cdef SIZE_t n_classes_node
    cdef SIZE_t n_classes_left
    cdef SIZE_t n_classes_right

    cdef DTYPE_t entropy_node
    cdef DTYPE_t entropy_left
    cdef DTYPE_t entropy_right

    cdef DTYPE_t frac
    cdef DTYPE_t gain = 0.0
    cdef DTYPE_t delta = 0.0

    cdef SIZE_t output
    cdef SIZE_t label

    cdef DTYPE_t count_label
    cdef SIZE_t label_index

    for output in range(n_outputs):
        # Reset the statistics for the output
        n_classes_node = 0
        n_classes_left = 0
        n_classes_right = 0

        entropy_node = 0.0
        entropy_left = 0.0
        entropy_right = 0.0

        for label in range(n_classes[output]):
            label_index = output * sum_stride + label

            # Compute the entropy and number of classes for the node
            count_label = sum_total[label_index] / n_samples_node

            if count_label > 0.0:
                n_classes_node += 1
                entropy_node -= count_label * log2(count_label)

            # Compute the entropy and number of classes for the left part
            count_label = sum_left[label_index] / n_samples_left

            if count_label > 0.0:
                n_classes_left += 1
                entropy_left -= count_label * log2(count_label)

            # Compute the entropy for the right part
            count_label = sum_right[label_index] / n_samples_right

            if count_label > 0.0:
                n_classes_right += 1
                entropy_right -= count_label * log2(count_label)

        gain += (entropy_node
                 - n_samples_left / n_samples_node * entropy_left
                 - n_samples_right / n_samples_node * entropy_right)

        delta += (log2(pow(3, n_classes_node) - 2)
                  - (n_classes_node * entropy_node
                     - n_classes_left * entropy_left
                     - n_classes_right * entropy_right))

    gain /= n_outputs
    delta /= n_outputs

    return gain > (log2(n_samples_node - 1) + delta) / n_samples_node


cdef DTYPE_t get_threshold(const DTYPE_t[:, :] X,
                           SIZE_t feature,
                           const SIZE_t[:, :] X_idx_sorted,
                           SIZE_t cut_point) nogil:
    """Get the threshold corresponding to a cut point."""
    return (X[X_idx_sorted[cut_point - 1, feature], feature]
            + X[X_idx_sorted[cut_point, feature], feature]) / 2


# =============================================================================
# Helpers
# =============================================================================

cdef void safe_calloc(PTR_t* pointer,
                      SIZE_t n_elements) nogil except *:
    """Safe allocate memory for a pointer."""
    cdef PTR_t temp = <PTR_t> calloc(n_elements, sizeof(pointer[0][0]))

    if temp == NULL:
        with gil:
            raise MemoryError("Could not allocate memory for "
                              "{0} elements.".format(n_elements))

    pointer[0] = temp


cdef void safe_realloc(PTR_t* pointer,
                       SIZE_t n_elements) nogil except *:
    """Safe re-allocate memory for a pointer."""
    cdef SIZE_t n_bytes = n_elements * sizeof(pointer[0][0])

    if n_bytes / sizeof(pointer[0][0]) != n_elements:
        # Overflow in the multiplication
        with gil:
            raise MemoryError("Could not allocate ({0} * {1}) bytes."
                              .format(n_elements, sizeof(pointer[0][0])))

    cdef PTR_t temp = <PTR_t> realloc(pointer[0], n_bytes)

    if temp == NULL:
        with gil:
            raise MemoryError("Could not allocate memory for "
                              "{0} bytes.".format(n_bytes))

    pointer[0] = temp


cdef inline DTYPE_t log2(DTYPE_t value) nogil:
    """Return the logarithm in base 2."""
    return log(value) / log(2)


# =============================================================================
# Stack data structure
# =============================================================================

cdef class Stack:
    """A stack data structure."""

    def __cinit__(self, SIZE_t capacity):
        self.capacity = capacity
        self.top = 0
        self.stack = NULL

        # Allocate memory for the stack
        safe_calloc(&self.stack, self.capacity)

    def __dealloc__(self):
        for top in range(self.top):
            # Deallocate all inner pointers
            free(self.stack[top].sum_total)

        free(self.stack)

    cdef int push(self,
                  SIZE_t start,
                  SIZE_t end,
                  DTYPE_t* sum_total,
                  SIZE_t sum_stride,
                  SIZE_t n_outputs) nogil except -1:
        """Push a new element onto the stack."""
        cdef SIZE_t n_elements = n_outputs * sum_stride
        cdef SIZE_t n_bytes = n_elements * sizeof(DTYPE_t)

        # Re-allocate memory for the stack
        if self.top >= self.capacity:
            self.capacity *= 2
            safe_realloc(&self.stack, self.capacity)

        self.stack[self.top].start = start
        self.stack[self.top].end = end

        safe_calloc(&self.stack[self.top].sum_total, n_elements)
        memcpy(self.stack[self.top].sum_total, sum_total, n_bytes)

        # Increment the stack pointer
        self.top += 1

        return 0

    cdef int pop(self,
                 StackRecord* stack_record,
                 SIZE_t sum_stride,
                 SIZE_t n_outputs) nogil except -1:
        """Remove the top element from the stack."""
        cdef SIZE_t n_elements = n_outputs * sum_stride
        cdef SIZE_t n_bytes = n_elements * sizeof(DTYPE_t)

        cdef DTYPE_t* sum_total = self.stack[self.top - 1].sum_total

        stack_record[0].start = self.stack[self.top - 1].start
        stack_record[0].end = self.stack[self.top - 1].end

        memcpy(stack_record[0].sum_total, sum_total, n_bytes)

        # The top record pointer can be safely deallocated
        free(sum_total)

        # Decrease the stack pointer
        self.top -= 1

        return 0
