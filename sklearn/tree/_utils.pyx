# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# Licence: BSD 3 clause

from libc.stdlib cimport free, malloc, realloc


# =============================================================================
# Stack data structure
# =============================================================================

cdef inline void copy_stack(StackRecord* a, StackRecord* b) nogil:
    """Assigns ``a := b`` for StackRecord. """
    a.start = b.start
    a.end = b.end
    a.depth = b.depth
    a.parent = b.parent
    a.is_left = b.is_left
    a.impurity = b.impurity


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

    cdef void push(self, SIZE_t start, SIZE_t end, SIZE_t depth, SIZE_t parent,
                   bint is_left, double impurity) nogil:
        """Push a new element onto the stack. """
        cdef SIZE_t top = self.top
        cdef StackRecord* stack = NULL

        # Resize if capacity not sufficient
        if top >= self.capacity:
            self.capacity *= 2
            self.stack_ = <StackRecord*> realloc(self.stack_, self.capacity * sizeof(StackRecord))

        stack = self.stack_
        stack[top].start = start
        stack[top].end = end
        stack[top].depth = depth
        stack[top].parent = parent
        stack[top].is_left = is_left
        stack[top].impurity = impurity

        # Increment stack pointer
        self.top = top + 1

    cdef int pop(self, StackRecord* res) nogil:
        """Remove the top element from the stack. """
        cdef SIZE_t top = self.top
        cdef StackRecord* stack = self.stack_

        if top <= 0:
            return 0

        copy_stack(res, stack + top - 1)
        self.top = top - 1

        return 1


# =============================================================================
# PriorityHeap data structure
# =============================================================================

cdef inline void copy_heap(PriorityHeapRecord* a, PriorityHeapRecord* b) nogil:
    """Assigns ``a := b``. """
    a.node_id = b.node_id
    a.start = b.start
    a.end = b.end
    a.pos = b.pos
    a.depth = b.depth
    a.is_leaf = b.is_leaf
    a.impurity = b.impurity
    a.improvement = b.improvement


cdef void swap_heap(PriorityHeapRecord* heap, SIZE_t a, SIZE_t b) nogil:
    """Swap record ``a`` and ``b`` in ``heap``. """
    cdef PriorityHeapRecord tmp
    copy_heap(&tmp, heap + a)
    copy_heap(heap + a, heap + b)
    copy_heap(heap + b, &tmp)


cdef void heapify_up(PriorityHeapRecord* heap, SIZE_t pos) nogil:
    """Restore heap invariant parent.improvement > child.improvement from ``pos`` upwards. """
    if pos == 0:
        return

    cdef SIZE_t parent_pos = (pos - 1) / 2

    if heap[parent_pos].improvement < heap[pos].improvement:
        swap_heap(heap, parent_pos, pos)
        heapify_up(heap, parent_pos)


cdef void heapify_down(PriorityHeapRecord* heap, SIZE_t pos, SIZE_t heap_length) nogil:
    """Restore heap invariant parent.improvement > children.improvement from ``pos`` downwards. """
    cdef SIZE_t left_pos = 2 * (pos + 1) - 1
    cdef SIZE_t right_pos = 2 * (pos + 1)
    cdef SIZE_t largest = pos

    if left_pos < heap_length and heap[left_pos].improvement > heap[largest].improvement:
        largest = left_pos
    if right_pos < heap_length and heap[right_pos].improvement > heap[largest].improvement:
        largest = right_pos

    if largest != pos:
        swap_heap(heap, pos, largest)
        heapify_down(heap, largest, heap_length)


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

    def __dealloc__(self):
        free(self.heap_)

    cdef bint is_empty(self) nogil:
        return self.heap_ptr <= 0

    cdef void push(self, SIZE_t node_id, SIZE_t start, SIZE_t end, SIZE_t pos,
                   SIZE_t depth, bint is_leaf, double improvement,
                   double impurity) nogil:
        """Push record on the priority heap. """
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef PriorityHeapRecord* heap = NULL

        # Resize if capacity not sufficient
        if heap_ptr >= self.capacity:
            self.capacity *= 2
            self.heap_ = <PriorityHeapRecord*> realloc(self.heap_,
                                                       self.capacity *
                                                       sizeof(PriorityHeapRecord))

        # Put element as last element of heap
        heap = self.heap_
        heap[heap_ptr].node_id = node_id
        heap[heap_ptr].start = start
        heap[heap_ptr].end = end
        heap[heap_ptr].pos = pos
        heap[heap_ptr].depth = depth
        heap[heap_ptr].is_leaf = is_leaf
        heap[heap_ptr].impurity = impurity
        heap[heap_ptr].improvement = improvement

        # Heapify up
        heapify_up(heap, heap_ptr)

        # Increase element count
        self.heap_ptr = heap_ptr + 1

    cdef int pop(self, PriorityHeapRecord* res) nogil:
        """Remove max element from the heap. """
        cdef SIZE_t heap_ptr = self.heap_ptr
        cdef PriorityHeapRecord* heap = self.heap_

        if heap_ptr <= 0:
            return 0

        # Take first element
        copy_heap(res, heap)

        # Put last element to the front
        swap_heap(heap, 0, heap_ptr - 1)

        # Restore heap invariant
        if heap_ptr > 1:
            heapify_down(heap, 0, heap_ptr - 1)

        self.heap_ptr = heap_ptr - 1

        return 1
