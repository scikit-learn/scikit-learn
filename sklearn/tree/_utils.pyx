# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#
#
# License: BSD 3 clause

from libc.stdlib cimport free
from libc.stdlib cimport malloc
from libc.stdlib cimport realloc
from libc.math cimport log as ln

import numpy as np
cimport numpy as np
np.import_array()

# =============================================================================
# Helper functions
# =============================================================================

cdef realloc_ptr safe_realloc(realloc_ptr* p, size_t nelems) nogil except *:
    # sizeof(realloc_ptr[0]) would be more like idiomatic C, but causes Cython
    # 0.20.1 to crash.
    cdef size_t nbytes = nelems * sizeof(p[0][0])
    if nbytes / sizeof(p[0][0]) != nelems:
        # Overflow in the multiplication
        with gil:
            raise MemoryError("could not allocate (%d * %d) bytes"
                              % (nelems, sizeof(p[0][0])))
    cdef realloc_ptr tmp = <realloc_ptr>realloc(p[0], nbytes)
    if tmp == NULL:
        with gil:
            raise MemoryError("could not allocate %d bytes" % nbytes)

    p[0] = tmp
    return tmp  # for convenience


def _realloc_test():
    # Helper for tests. Tries to allocate <size_t>(-1) / 2 * sizeof(size_t)
    # bytes, which will always overflow.
    cdef SIZE_t* p = NULL
    safe_realloc(&p, <size_t>(-1) / 2)
    if p != NULL:
        free(p)
        assert False


# rand_r replacement using a 32bit XorShift generator
# See https://www.jstatsoft.org/v08/i14/paper for details
cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)

    return seed[0] % (<UINT32_t>RAND_R_MAX + 1)


cdef inline np.ndarray sizet_ptr_to_ndarray(SIZE_t* data, SIZE_t size):
    """Return copied data as 1D numpy array of intp's."""
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> size
    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INTP, data).copy()


cdef inline SIZE_t rand_int(SIZE_t low, SIZE_t high,
                            UINT32_t* random_state) nogil:
    """Generate a random integer in [low; end)."""
    return low + our_rand_r(random_state) % (high - low)


cdef inline double rand_uniform(double low, double high,
                                UINT32_t* random_state) nogil:
    """Generate a random double in [low; high)."""
    return ((high - low) * <double> our_rand_r(random_state) /
            <double> RAND_R_MAX) + low


cdef inline double log(double x) nogil:
    return ln(x) / ln(2.0)


# =============================================================================
# Stack data structure
# =============================================================================

cdef class Stack:
    """A LIFO data structure.

    Attributes
    ----------
    capacity : SIZE_t
        The elements the stack can hold; if more added then ``self.stack_``
        needs to be resized.

    top : SIZE_t
        The number of elements currently on the stack.

    stack : StackRecord pointer
        The stack of records (upward in the stack corresponds to the right).
    """

    def __cinit__(self, SIZE_t capacity):
        self.capacity = capacity
        self.top = 0
        self.stack_ = <StackRecord*> malloc(capacity * sizeof(StackRecord))

    def __dealloc__(self):
        free(self.stack_)

    cdef bint is_empty(self) nogil:
        return self.top <= 0

    cdef int push(self, SIZE_t start, SIZE_t end, SIZE_t depth, SIZE_t parent,
                  bint is_left, double impurity,
                  SIZE_t n_constant_features) nogil except -1:
        """Push a new element onto the stack.

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef SIZE_t top = self.top
        cdef StackRecord* stack = NULL

        # Resize if capacity not sufficient
        if top >= self.capacity:
            self.capacity *= 2
            # Since safe_realloc can raise MemoryError, use `except -1`
            safe_realloc(&self.stack_, self.capacity)

        stack = self.stack_
        stack[top].start = start
        stack[top].end = end
        stack[top].depth = depth
        stack[top].parent = parent
        stack[top].is_left = is_left
        stack[top].impurity = impurity
        stack[top].n_constant_features = n_constant_features

        # Increment stack pointer
        self.top = top + 1
        return 0

    cdef int pop(self, StackRecord* res) nogil:
        """Remove the top element from the stack and copy to ``res``.

        Returns 0 if pop was successful (and ``res`` is set); -1
        otherwise.
        """
        cdef SIZE_t top = self.top
        cdef StackRecord* stack = self.stack_

        if top <= 0:
            return -1

        res[0] = stack[top - 1]
        self.top = top - 1

        return 0


# =============================================================================
# PriorityHeap data structure
# =============================================================================

cdef class PriorityHeap:
    """A priority queue implemented as a binary heap.

    The heap invariant is that the impurity improvement of the parent record
    is larger then the impurity improvement of the children.

    Attributes
    ----------
    capacity : SIZE_t
        The capacity of the heap

    heap_ptr : SIZE_t
        The water mark of the heap; the heap grows from left to right in the
        array ``heap_``. The following invariant holds ``heap_ptr < capacity``.

    heap_ : PriorityHeapRecord*
        The array of heap records. The maximum element is on the left;
        the heap grows from left to right
    """

    def __cinit__(self, SIZE_t capacity):
        self.capacity = capacity
        self.heap_ptr = 0
        safe_realloc(&self.heap_, capacity)

    def __dealloc__(self):
        free(self.heap_)

    cdef bint is_empty(self) nogil:
        return self.heap_ptr <= 0

    cdef void heapify_up(self, PriorityHeapRecord* heap, SIZE_t pos) nogil:
        """Restore heap invariant parent.improvement > child.improvement from
           ``pos`` upwards. """
        if pos == 0:
            return

        cdef SIZE_t parent_pos = (pos - 1) / 2

        if heap[parent_pos].improvement < heap[pos].improvement:
            heap[parent_pos], heap[pos] = heap[pos], heap[parent_pos]
            self.heapify_up(heap, parent_pos)

    cdef void heapify_down(self, PriorityHeapRecord* heap, SIZE_t pos,
                           SIZE_t heap_length) nogil:
        """Restore heap invariant parent.improvement > children.improvement from
           ``pos`` downwards. """
        cdef SIZE_t left_pos = 2 * (pos + 1) - 1
        cdef SIZE_t right_pos = 2 * (pos + 1)
        cdef SIZE_t largest = pos

        if (left_pos < heap_length and
                heap[left_pos].improvement > heap[largest].improvement):
            largest = left_pos

        if (right_pos < heap_length and
                heap[right_pos].improvement > heap[largest].improvement):
            largest = right_pos

        if largest != pos:
            heap[pos], heap[largest] = heap[largest], heap[pos]
            self.heapify_down(heap, largest, heap_length)

    cdef int push(self, SIZE_t node_id, SIZE_t start, SIZE_t end, SIZE_t pos,
                  SIZE_t depth, bint is_leaf, double improvement,
                  double impurity, double impurity_left,
                  double impurity_right) nogil except -1:
        """Push record on the priority heap.

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef PriorityHeapRecord* heap = NULL

        # Resize if capacity not sufficient
        if heap_ptr >= self.capacity:
            self.capacity *= 2
            # Since safe_realloc can raise MemoryError, use `except -1`
            safe_realloc(&self.heap_, self.capacity)

        # Put element as last element of heap
        heap = self.heap_
        heap[heap_ptr].node_id = node_id
        heap[heap_ptr].start = start
        heap[heap_ptr].end = end
        heap[heap_ptr].pos = pos
        heap[heap_ptr].depth = depth
        heap[heap_ptr].is_leaf = is_leaf
        heap[heap_ptr].impurity = impurity
        heap[heap_ptr].impurity_left = impurity_left
        heap[heap_ptr].impurity_right = impurity_right
        heap[heap_ptr].improvement = improvement

        # Heapify up
        self.heapify_up(heap, heap_ptr)

        # Increase element count
        self.heap_ptr = heap_ptr + 1
        return 0

    cdef int pop(self, PriorityHeapRecord* res) nogil:
        """Remove max element from the heap. """
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef PriorityHeapRecord* heap = self.heap_

        if heap_ptr <= 0:
            return -1

        # Take first element
        res[0] = heap[0]

        # Put last element to the front
        heap[0], heap[heap_ptr - 1] = heap[heap_ptr - 1], heap[0]

        # Restore heap invariant
        if heap_ptr > 1:
            self.heapify_down(heap, 0, heap_ptr - 1)

        self.heap_ptr = heap_ptr - 1

        return 0

# =============================================================================
# WeightedPQueue data structure
# =============================================================================

cdef class WeightedPQueue:
    """A priority queue class, always sorted in increasing order.

    Attributes
    ----------
    capacity : SIZE_t
        The capacity of the priority queue.

    array_ptr : SIZE_t
        The water mark of the priority queue; the priority queue grows from
        left to right in the array ``array_``. ``array_ptr`` is always
        less than ``capacity``.

    array_ : WeightedPQueueRecord*
        The array of priority queue records. The minimum element is on the
        left at index 0, and the maximum element is on the right at index
        ``array_ptr-1``.
    """

    def __cinit__(self, SIZE_t capacity):
        self.capacity = capacity
        self.array_ptr = 0
        safe_realloc(&self.array_, capacity)

    def __dealloc__(self):
        free(self.array_)

    cdef int reset(self) nogil except -1:
        """Reset the WeightedPQueue to its state at construction

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        self.array_ptr = 0
        # Since safe_realloc can raise MemoryError, use `except *`
        safe_realloc(&self.array_, self.capacity)
        return 0

    cdef bint is_empty(self) nogil:
        return self.array_ptr <= 0

    cdef SIZE_t size(self) nogil:
        return self.array_ptr

    cdef int push(self, DOUBLE_t data, DOUBLE_t weight) nogil except -1:
        """Push record on the array.

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef SIZE_t array_ptr = self.array_ptr
        cdef WeightedPQueueRecord* array = NULL
        cdef SIZE_t i

        # Resize if capacity not sufficient
        if array_ptr >= self.capacity:
            self.capacity *= 2
            # Since safe_realloc can raise MemoryError, use `except -1`
            safe_realloc(&self.array_, self.capacity)

        # Put element as last element of array
        array = self.array_
        array[array_ptr].data = data
        array[array_ptr].weight = weight

        # bubble last element up according until it is sorted
        # in ascending order
        i = array_ptr
        while(i != 0 and array[i].data < array[i-1].data):
            array[i], array[i-1] = array[i-1], array[i]
            i -= 1

        # Increase element count
        self.array_ptr = array_ptr + 1
        return 0

    cdef int remove(self, DOUBLE_t data, DOUBLE_t weight) nogil:
        """Remove a specific value/weight record from the array.
        Returns 0 if successful, -1 if record not found."""
        cdef SIZE_t array_ptr = self.array_ptr
        cdef WeightedPQueueRecord* array = self.array_
        cdef SIZE_t idx_to_remove = -1
        cdef SIZE_t i

        if array_ptr <= 0:
            return -1

        # find element to remove
        for i in range(array_ptr):
            if array[i].data == data and array[i].weight == weight:
                idx_to_remove = i
                break

        if idx_to_remove == -1:
            return -1

        # shift the elements after the removed element
        # to the left.
        for i in range(idx_to_remove, array_ptr-1):
            array[i] = array[i+1]

        self.array_ptr = array_ptr - 1
        return 0

    cdef int pop(self, DOUBLE_t* data, DOUBLE_t* weight) nogil:
        """Remove the top (minimum) element from array.
        Returns 0 if successful, -1 if nothing to remove."""
        cdef SIZE_t array_ptr = self.array_ptr
        cdef WeightedPQueueRecord* array = self.array_
        cdef SIZE_t i

        if array_ptr <= 0:
            return -1

        data[0] = array[0].data
        weight[0] = array[0].weight

        # shift the elements after the removed element
        # to the left.
        for i in range(0, array_ptr-1):
            array[i] = array[i+1]

        self.array_ptr = array_ptr - 1
        return 0

    cdef int peek(self, DOUBLE_t* data, DOUBLE_t* weight) nogil:
        """Write the top element from array to a pointer.
        Returns 0 if successful, -1 if nothing to write."""
        cdef WeightedPQueueRecord* array = self.array_
        if self.array_ptr <= 0:
            return -1
        # Take first value
        data[0] = array[0].data
        weight[0] = array[0].weight
        return 0

    cdef DOUBLE_t get_weight_from_index(self, SIZE_t index) nogil:
        """Given an index between [0,self.current_capacity], access
        the appropriate heap and return the requested weight"""
        cdef WeightedPQueueRecord* array = self.array_

        # get weight at index
        return array[index].weight

    cdef DOUBLE_t get_value_from_index(self, SIZE_t index) nogil:
        """Given an index between [0,self.current_capacity], access
        the appropriate heap and return the requested value"""
        cdef WeightedPQueueRecord* array = self.array_

        # get value at index
        return array[index].data

# =============================================================================
# WeightedMedianCalculator data structure
# =============================================================================

cdef class WeightedMedianCalculator:
    """A class to handle calculation of the weighted median from streams of
    data. To do so, it maintains a parameter ``k`` such that the sum of the
    weights in the range [0,k) is greater than or equal to half of the total
    weight. By minimizing the value of ``k`` that fulfills this constraint,
    calculating the median is done by either taking the value of the sample
    at index ``k-1`` of ``samples`` (samples[k-1].data) or the average of
    the samples at index ``k-1`` and ``k`` of ``samples``
    ((samples[k-1] + samples[k]) / 2).

    Attributes
    ----------
    initial_capacity : SIZE_t
        The initial capacity of the WeightedMedianCalculator.

    samples : WeightedPQueue
        Holds the samples (consisting of values and their weights) used in the
        weighted median calculation.

    total_weight : DOUBLE_t
        The sum of the weights of items in ``samples``. Represents the total
        weight of all samples used in the median calculation.

    k : SIZE_t
        Index used to calculate the median.

    sum_w_0_k : DOUBLE_t
        The sum of the weights from samples[0:k]. Used in the weighted
        median calculation; minimizing the value of ``k`` such that
        ``sum_w_0_k`` >= ``total_weight / 2`` provides a mechanism for
        calculating the median in constant time.

    """

    def __cinit__(self, SIZE_t initial_capacity):
        self.initial_capacity = initial_capacity
        self.samples = WeightedPQueue(initial_capacity)
        self.total_weight = 0
        self.k = 0
        self.sum_w_0_k = 0

    cdef SIZE_t size(self) nogil:
        """Return the number of samples in the
        WeightedMedianCalculator"""
        return self.samples.size()

    cdef int reset(self) nogil except -1:
        """Reset the WeightedMedianCalculator to its state at construction

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        # samples.reset (WeightedPQueue.reset) uses safe_realloc, hence
        # except -1
        self.samples.reset()
        self.total_weight = 0
        self.k = 0
        self.sum_w_0_k = 0
        return 0

    cdef int push(self, DOUBLE_t data, DOUBLE_t weight) nogil except -1:
        """Push a value and its associated weight to the WeightedMedianCalculator

        Return -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        cdef int return_value
        cdef DOUBLE_t original_median

        if self.size() != 0:
            original_median = self.get_median()
        # samples.push (WeightedPQueue.push) uses safe_realloc, hence except -1
        return_value = self.samples.push(data, weight)
        self.update_median_parameters_post_push(data, weight,
                                                original_median)
        return return_value

    cdef int update_median_parameters_post_push(
            self, DOUBLE_t data, DOUBLE_t weight,
            DOUBLE_t original_median) nogil:
        """Update the parameters used in the median calculation,
        namely `k` and `sum_w_0_k` after an insertion"""

        # trivial case of one element.
        if self.size() == 1:
            self.k = 1
            self.total_weight = weight
            self.sum_w_0_k = self.total_weight
            return 0

        # get the original weighted median
        self.total_weight += weight

        if data < original_median:
            # inserting below the median, so increment k and
            # then update self.sum_w_0_k accordingly by adding
            # the weight that was added.
            self.k += 1
            # update sum_w_0_k by adding the weight added
            self.sum_w_0_k += weight

            # minimize k such that sum(W[0:k]) >= total_weight / 2
            # minimum value of k is 1
            while(self.k > 1 and ((self.sum_w_0_k -
                                   self.samples.get_weight_from_index(self.k-1))
                                  >= self.total_weight / 2.0)):
                self.k -= 1
                self.sum_w_0_k -= self.samples.get_weight_from_index(self.k)
            return 0

        if data >= original_median:
            # inserting above or at the median
            # minimize k such that sum(W[0:k]) >= total_weight / 2
            while(self.k < self.samples.size() and
                  (self.sum_w_0_k < self.total_weight / 2.0)):
                self.k += 1
                self.sum_w_0_k += self.samples.get_weight_from_index(self.k-1)
            return 0

    cdef int remove(self, DOUBLE_t data, DOUBLE_t weight) nogil:
        """Remove a value from the MedianHeap, removing it
        from consideration in the median calculation
        """
        cdef int return_value
        cdef DOUBLE_t original_median

        if self.size() != 0:
            original_median = self.get_median()

        return_value = self.samples.remove(data, weight)
        self.update_median_parameters_post_remove(data, weight,
                                                  original_median)
        return return_value

    cdef int pop(self, DOUBLE_t* data, DOUBLE_t* weight) nogil:
        """Pop a value from the MedianHeap, starting from the
        left and moving to the right.
        """
        cdef int return_value
        cdef double original_median

        if self.size() != 0:
            original_median = self.get_median()

        # no elements to pop
        if self.samples.size() == 0:
            return -1

        return_value = self.samples.pop(data, weight)
        self.update_median_parameters_post_remove(data[0],
                                                  weight[0],
                                                  original_median)
        return return_value

    cdef int update_median_parameters_post_remove(
            self, DOUBLE_t data, DOUBLE_t weight,
            double original_median) nogil:
        """Update the parameters used in the median calculation,
        namely `k` and `sum_w_0_k` after a removal"""
        # reset parameters because it there are no elements
        if self.samples.size() == 0:
            self.k = 0
            self.total_weight = 0
            self.sum_w_0_k = 0
            return 0

        # trivial case of one element.
        if self.samples.size() == 1:
            self.k = 1
            self.total_weight -= weight
            self.sum_w_0_k = self.total_weight
            return 0

        # get the current weighted median
        self.total_weight -= weight

        if data < original_median:
            # removing below the median, so decrement k and
            # then update self.sum_w_0_k accordingly by subtracting
            # the removed weight

            self.k -= 1
            # update sum_w_0_k by removing the weight at index k
            self.sum_w_0_k -= weight

            # minimize k such that sum(W[0:k]) >= total_weight / 2
            # by incrementing k and updating sum_w_0_k accordingly
            # until the condition is met.
            while(self.k < self.samples.size() and
                  (self.sum_w_0_k < self.total_weight / 2.0)):
                self.k += 1
                self.sum_w_0_k += self.samples.get_weight_from_index(self.k-1)
            return 0

        if data >= original_median:
            # removing above the median
            # minimize k such that sum(W[0:k]) >= total_weight / 2
            while(self.k > 1 and ((self.sum_w_0_k -
                                   self.samples.get_weight_from_index(self.k-1))
                                  >= self.total_weight / 2.0)):
                self.k -= 1
                self.sum_w_0_k -= self.samples.get_weight_from_index(self.k)
            return 0

    cdef DOUBLE_t get_median(self) nogil:
        """Write the median to a pointer, taking into account
        sample weights."""
        if self.sum_w_0_k == (self.total_weight / 2.0):
            # split median
            return (self.samples.get_value_from_index(self.k) +
                    self.samples.get_value_from_index(self.k-1)) / 2.0
        if self.sum_w_0_k > (self.total_weight / 2.0):
            # whole median
            return self.samples.get_value_from_index(self.k-1)
