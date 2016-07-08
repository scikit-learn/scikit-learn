# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

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
from libc.stdlib cimport calloc
from libc.stdlib cimport realloc
from libc.math cimport log as ln

import numpy as np
cimport numpy as np
np.import_array()

# =============================================================================
# Helper functions
# =============================================================================

cdef realloc_ptr safe_realloc(realloc_ptr* p, size_t nelems) except *:
    # sizeof(realloc_ptr[0]) would be more like idiomatic C, but causes Cython
    # 0.20.1 to crash.
    cdef size_t nbytes = nelems * sizeof(p[0][0])
    if nbytes / sizeof(p[0][0]) != nelems:
        # Overflow in the multiplication
        raise MemoryError("could not allocate (%d * %d) bytes"
                          % (nelems, sizeof(p[0][0])))
    cdef realloc_ptr tmp = <realloc_ptr>realloc(p[0], nbytes)
    if tmp == NULL:
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
# See http://www.jstatsoft.org/v08/i14/paper for details
cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)

    return seed[0] % (<UINT32_t>RAND_R_MAX + 1)


cdef inline np.ndarray sizet_ptr_to_ndarray(SIZE_t* data, SIZE_t size):
    """Encapsulate data into a 1D numpy array of intp's."""
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> size
    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INTP, data)


cdef inline SIZE_t rand_int(SIZE_t low, SIZE_t high,
                            UINT32_t* random_state) nogil:
    """Generate a random integer in [0; end)."""
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
        if self.stack_ == NULL:
            raise MemoryError()

    def __dealloc__(self):
        free(self.stack_)

    cdef bint is_empty(self) nogil:
        return self.top <= 0

    cdef int push(self, SIZE_t start, SIZE_t end, SIZE_t depth, SIZE_t parent,
                  bint is_left, double impurity,
                  SIZE_t n_constant_features) nogil:
        """Push a new element onto the stack.

        Returns 0 if successful; -1 on out of memory error.
        """
        cdef SIZE_t top = self.top
        cdef StackRecord* stack = NULL

        # Resize if capacity not sufficient
        if top >= self.capacity:
            self.capacity *= 2
            stack = <StackRecord*> realloc(self.stack_,
                                           self.capacity * sizeof(StackRecord))
            if stack == NULL:
                # no free; __dealloc__ handles that
                return -1
            self.stack_ = stack

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
        self.heap_ = <PriorityHeapRecord*> malloc(capacity * sizeof(PriorityHeapRecord))
        if self.heap_ == NULL:
            raise MemoryError()

    def __dealloc__(self):
        free(self.heap_)

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

    cdef bint is_empty(self) nogil:
        return self.heap_ptr <= 0

    cdef int push(self, SIZE_t node_id, SIZE_t start, SIZE_t end, SIZE_t pos,
                  SIZE_t depth, bint is_leaf, double improvement,
                  double impurity, double impurity_left,
                  double impurity_right) nogil:
        """Push record on the priority heap.

        Returns 0 if successful; -1 on out of memory error.
        """
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef PriorityHeapRecord* heap = NULL

        # Resize if capacity not sufficient
        if heap_ptr >= self.capacity:
            self.capacity *= 2
            heap = <PriorityHeapRecord*> realloc(self.heap_,
                                                 self.capacity *
                                                 sizeof(PriorityHeapRecord))
            if heap == NULL:
                # no free; __dealloc__ handles that
                return -1
            self.heap_ = heap

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
# MinMaxHeap data structure
# =============================================================================

cdef class MinMaxHeap:
    """A priority queue implemented as a binary heap.

    The heap invariant is that the impurity improvement of the parent record is
    larger then the impurity improvement of the children. The MinHeap is
    essentially an array sorted in ascending order, and a MaxHeap is an array
    sorted in descending order.

    Attributes
    ----------
    capacity : SIZE_t
        The capacity of the heap

    heap_ptr : SIZE_t
        The water mark of the heap; the heap grows from left to right in the
        array ``heap_``. heap_ptr is always less than capacity.

    heap_ : MinMaxHeapRecord*
        The array of heap records. The maximum element is on the left;
        the heap grows from left to right

    mode : bint
        The mode of the heap. When the value of the ``mode`` parameter passed
        in at construction is ``max``, the heap is a Max-Heap and mode is set
        to 1. When the value of the ``mode`` parameter passed in at
        construction is not ``max``, the heap is a Min-Heap and mode is set
        to 0.
    """

    def __cinit__(self, SIZE_t capacity, str mode):
        self.capacity = capacity
        if mode == "max":
            self.mode = 1
        else:
            self.mode = 0

        self.heap_ptr = 0

        self.heap_ = <MinMaxHeapRecord*> calloc(capacity, sizeof(MinMaxHeapRecord))
        if self.heap_ == NULL:
            raise MemoryError()

    def __dealloc__(self):
        free(self.heap_)

    cdef bint is_empty(self) nogil:
        return self.heap_ptr <= 0

    cdef SIZE_t size(self) nogil:
        return self.heap_ptr

    cdef int push(self, DOUBLE_t data) nogil:
        """Push record on the priority heap.

        Returns 0 if successful; -1 on out of memory error.
        """
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef SIZE_t i
        cdef MinMaxHeapRecord* heap = NULL

        # Resize if capacity not sufficient
        if heap_ptr >= self.capacity:
            self.capacity *= 2
            heap = <MinMaxHeapRecord*> realloc(self.heap_,
                                               self.capacity *
                                               sizeof(MinMaxHeapRecord))
            if heap == NULL:
                # no free; __dealloc__ handles that
                return -1
            self.heap_ = heap

        # Put element as last element of heap
        heap = self.heap_
        heap[heap_ptr].data = data

        # bubble last element up according to mode
        # max heap, sorted in descending order
        i = heap_ptr
        if self.mode == 1:
            while(i != 0 and heap[i].data > heap[i-1].data):
                heap[i], heap[i-1] = heap[i-1], heap[i]
                i = i-1

        # min heap, sorted in ascending order
        else:
            while(i != 0 and heap[i].data < heap[i-1].data):
                heap[i], heap[i-1] = heap[i-1], heap[i]
                i = i-1

        # Increase element count
        self.heap_ptr = heap_ptr + 1
        return 0

    cdef int remove(self, DOUBLE_t value) nogil:
        """Remove a specific value from heap"""
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef MinMaxHeapRecord* heap = self.heap_
        cdef SIZE_t idx_to_remove = -1
        cdef SIZE_t i

        if heap_ptr <= 0:
            return -1

        # find element to remove
        for i in range(heap_ptr):
            if heap[i].data == value:
                idx_to_remove = i
                break
        # should we throw an error if the element isn't found?
        # it shouldn't happen, but better to fail noisily...?

        # move after the removed element over by one
        for i in range(idx_to_remove, heap_ptr-1):
            heap[i] = heap[i+1]

        self.heap_ptr = heap_ptr - 1
        return 0

    cdef int pop(self, DOUBLE_t* res) nogil:
        """Remove top element from heap."""
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef MinMaxHeapRecord* heap = self.heap_
        cdef SIZE_t i

        if heap_ptr <= 0:
            return -1

        # Take first element
        res[0] = heap[0].data

        # move after the removed element over by one
        for i in range(0, heap_ptr-1):
            heap[i] = heap[i+1]

        self.heap_ptr = heap_ptr - 1
        return 0

    cdef int peek(self, DOUBLE_t* res) nogil:
        """Write the top element from heap to a pointer."""
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef MinMaxHeapRecord* heap = self.heap_
        if heap_ptr <= 0:
            return -1
        # Take first value
        res[0] = heap[0].data
        return 0

# =============================================================================
# MedianHeap data structure
# =============================================================================

cdef class MedianHeap:

    def __cinit__(self, SIZE_t initial_capacity):
        self.initial_capacity = initial_capacity
        self.current_capacity = 0
        self.left_max_heap = MinMaxHeap(initial_capacity, "max")
        self.right_min_heap = MinMaxHeap(initial_capacity, "min")

    cdef SIZE_t size(self) nogil:
        return self.current_capacity

    cdef int push(self, DOUBLE_t data) nogil:
        """Push a value to the MedianHeap to be considered
        in the median calculation
        """
        cdef double current_median
        cdef int return_value

        if self.current_capacity == 0:
            return_value = self.left_max_heap.push(data)
        else:
            self.get_median(&current_median)
            if current_median <= data:
                # data is greater than or equal to current median, so it goes on min heap
                return_value = self.right_min_heap.push(data)
            else:
                # data is less than current median, so it goes on max heap
                return_value = self.left_max_heap.push(data)
        self.rebalance()
        self.current_capacity += 1
        return return_value

    cdef int remove(self, DOUBLE_t data) nogil:
        """Remove a value from the MedianHeap, removing it
        from consideration in the median calculation
        """
        cdef double current_median
        cdef int return_value

        self.get_median(&current_median)
        if current_median == data:
            # data is the same value as current median, it is in
            # the bigger one
            if self.right_min_heap.size() > self.left_max_heap.size():
                # it is in the right
                return_value = self.right_min_heap.remove(data)
            else:
                # it is in the left
                return_value = self.left_max_heap.remove(data)
        elif current_median < data:
            # data is greater than or equal to current median, so it is on min heap
            return_value = self.right_min_heap.remove(data)
        else:
            # data is less than current median, so it is on max heap
            return_value = self.left_max_heap.remove(data)
        self.rebalance()
        self.current_capacity -= 1
        return return_value

    cdef int pop(self, DOUBLE_t* res) nogil:
        """Pop a value from the MedianHeap, starting from the
        left and moving to the right.
        """
        cdef int return_value

        # no elements to pop
        if self.current_capacity == 0:
            return -1

        if self.left_max_heap.size() != 0:
            # pop from the left
            return_value = self.left_max_heap.pop(res)
        elif self.right_min_heap.size() != 0:
            # pop from right
            return_value = self.right_min_heap.pop(res)
        else:
            return -1
        self.rebalance()
        self.current_capacity -= 1
        return return_value

    cdef int get_median(self, double* data) nogil:
        """Return the current median"""
        if self.current_capacity == 0:
            return -1

        cdef SIZE_t left_max_heap_size = self.left_max_heap.size()
        cdef SIZE_t right_min_heap_size = self.right_min_heap.size()
        cdef DOUBLE_t left_max_heap_median
        cdef DOUBLE_t right_min_heap_median

        if self.current_capacity < 2:
            # there is only one thing, so set the median to be that
            if left_max_heap_size >= 1:
                self.left_max_heap.peek(&left_max_heap_median)
                data[0] = left_max_heap_median
            else:
                self.right_min_heap.peek(&right_min_heap_median)
                data[0] = right_min_heap_median
            return 0
        self.left_max_heap.peek(&left_max_heap_median)
        self.right_min_heap.peek(&right_min_heap_median)

        if left_max_heap_size == right_min_heap_size:
            # take the average of the two
            data[0] = (left_max_heap_median +
                            right_min_heap_median) / 2.0
        elif left_max_heap_size > right_min_heap_size:
            # left max heap larger, so median is at its' top
            data[0] = left_max_heap_median
        else:
            # right min heap is larger, so median is at its' top
            data[0] = right_min_heap_median
        return 0

    cdef int rebalance(self) nogil:
        """Rebalance the left max heap and the left min heap to have a
        one element or less difference in size"""
        cdef SIZE_t left_max_heap_size = self.left_max_heap.size()
        cdef SIZE_t right_min_heap_size = self.right_min_heap.size()
        cdef SIZE_t size_difference = left_max_heap_size - right_min_heap_size
        cdef DOUBLE_t popped
        cdef SIZE_t i

        if size_difference >= -1 and size_difference <= 1:
            # no balancing needed
            return 0

        if size_difference > 1:
            # left max heap bigger
            for i in range(0, size_difference - 1):
                # pop from left max heap and push into right min heap
                self.left_max_heap.pop(&popped)
                self.right_min_heap.push(popped)
        else:
            # right min heap bigger
            for i in range(0, (size_difference * -1) - 1):
                # pop from right min heap and push into left max heap
                self.right_min_heap.pop(&popped)
                self.left_max_heap.push(popped)
        return 0
