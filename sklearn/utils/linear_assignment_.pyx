# coding: UTF-8
"""
    @author: samuel ko
    @target: Transport linear_assignment on Cython.
    @notice: There is still many places to optimize, could anyone give me some advise to
             make this cython version of linear assignment more quickly and useful?

             Appreciate it!
"""

# Author: Gao Daiheng
# License: BSD
#
# cython: language_level=3

from __future__ import division
cimport cython
cimport numpy as np
import numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
# np.import_array()


DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
TTYPE = np.int64
ctypedef np.int64_t TTYPE_t


def linear_assignment(X):
    indices = _hungarian(X).tolist()
    indices.sort()
    # Re-force dtype to ints in case of empty list
    indices = np.array(indices, dtype=int)
    # Make sure the array is 2D with 2 columns.
    # This is needed when dealing with an empty list
    indices.shape = (-1, 2)
    return indices


def _hungarian(cost_matrix):
    state = _HungarianState(np.atleast_2d(cost_matrix))

    # No need to bother with assignments if one of the dimensions
    # of the cost matrix is zero-length.
    step = None if 0 in cost_matrix.shape else _step1

    while step is not None:
        step = step(state)

    # Look for the starred columns
    results = np.array(np.where(state.marked.base == int(1))).T

    # We need to swap the columns because we originally
    # did a transpose on the input cost matrix.
    if state.transposed:
        results = results[:, ::-1]

    return results


cdef class _HungarianState(object):
    """State of one execution of the Hungarian algorithm.

    Parameters
    ----------
    cost_matrix : 2D matrix
        The cost matrix. Does not have to be square.
    """

    cdef public double[:, :] C
    cdef public int transposed
    cdef int Z0_r, Z0_c
    cdef int n, m
    cdef public int [:] row_uncovered, col_uncovered
    cdef public int[:, :] path, marked


    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, np.ndarray[double, ndim=2] cost_matrix):

        cdef np.npy_bool transposed = (cost_matrix.shape[1] < cost_matrix.shape[0])
        if transposed:
            self.C = (cost_matrix.T).copy()
        else:
            self.C = cost_matrix.copy()

        self.transposed = transposed
        self.C = np.asarray(self.C, dtype=np.double)
        # At this point, m >= n.
        self.n = self.C.shape[0]
        self.m = self.C.shape[1]
        self.row_uncovered = np.ones(self.n, dtype=np.int32)
        self.col_uncovered = np.ones(self.m, dtype=np.int32)


        self.Z0_r = 0
        self.Z0_c = 0
        self.path = np.zeros((self.n + self.m, 2), dtype=np.int32)
        self.marked = np.zeros((self.n, self.m), dtype=np.int32)

    def _clear_covers(self):
        """Clear all covered matrix cells"""
        self.row_uncovered[:] = 1
        self.col_uncovered[:] = 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef _step1(_HungarianState state):
    """Steps 1 and 2 in the Wikipedia page."""

    # Step1: For each row of the matrix, find the smallest element and
    # subtract it from every element in its row.
    state.C -= state.C.base.min(axis=1)[:, np.newaxis]
    # Step2: Find a zero (Z) in the resulting matrix. If there is no
    # starred zero in its row or column, star Z. Repeat for each element
    # in the matrix.
    for i, j in zip(*np.where(state.C == np.float64(0))):
        if state.col_uncovered[j] and state.row_uncovered[i]:
            state.marked[i, j] = 1
            state.col_uncovered[j] = 0
            state.row_uncovered[i] = 0

    state._clear_covers()
    return _step3

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _step3(_HungarianState state):
    """
    Cover each column containing a starred zero. If n columns are covered,
    the starred zeros describe a complete set of unique assignments.
    In this case, Go to DONE, otherwise, Go to Step 4.
    """
    cdef np.ndarray[int, ndim=2] marked = (state.marked.base == int(1)).astype(dtype=np.int32)
    cdef np.ndarray[int, ndim=1] tmp_arr = np.any(marked, axis=0).astype(np.int32)

    for idx, i in enumerate(tmp_arr):
        if i == 1:
            state.col_uncovered[idx] = 0

    if marked.sum() < state.n:
        return _step4

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _step4(_HungarianState state):
    """
    Find a noncovered zero and prime it. If there is no starred zero
    in the row containing this primed zero, Go to Step 5. Otherwise,
    cover this row and uncover the column containing the starred
    zero. Continue in this manner until there are no uncovered zeros
    left. Save the smallest uncovered value and Go to Step 6.
    """
    # We convert to int as numpy operations are faster on int
    cdef np.ndarray[int, ndim=2] C = (state.C.base == np.double(0)).astype(np.int32)
    cdef np.ndarray[int, ndim=2] covered_C = C * state.row_uncovered.base[:, np.newaxis]
    covered_C *= state.col_uncovered.base

    cdef int row, col, star_col
    while True:
        # Find an uncovered zero
        row, col = np.unravel_index(np.argmax(covered_C), (state.n, state.m))
        if covered_C[row, col] == int(0):
            return _step6
        else:
            state.marked[row, col] = 2
            # Find the first starred element in the row
            star_col = np.argmax(state.marked.base[row, :] == int(1))

            if state.marked.base[row, star_col] != int(1):
                # Could not find one
                state.Z0_r = row
                state.Z0_c = col
                return _step5
            else:
                col = star_col
                state.row_uncovered[row] = 0
                state.col_uncovered[col] = 1
                covered_C[:, col] = C[:, col] * state.row_uncovered.base
                covered_C[row, :] = 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _step5(_HungarianState state):
    """
    Construct a series of alternating primed and starred zeros as follows.
    Let Z0 represent the uncovered primed zero found in Step 4.
    Let Z1 denote the starred zero in the column of Z0 (if any).
    Let Z2 denote the primed zero in the row of Z1 (there will always be one).
    Continue until the series terminates at a primed zero that has no starred
    zero in its column. Unstar each starred zero of the series, star each
    primed zero of the series, erase all primes and uncover every line in the
    matrix. Return to Step 3
    """
    cdef int count = 0
    state.path[count, 0] = state.Z0_r
    state.path[count, 1] = state.Z0_c

    cdef int row, col
    while True:
        # Find the first starred element in the col defined by
        # the path.
        row = np.argmax(state.marked.base[:, state.path[count, 1]] == int(1))

        if state.marked.base[row, state.path[count, 1]] != int(1):
            # Could not find one
            break
        else:
            count += 1
            state.path[count, 0] = row
            state.path[count, 1] = state.path[count - 1, 1]

        # Find the first prime element in the row defined by the
        # first path step
        col = np.argmax(state.marked.base[state.path[count, 0], :] == int(2))

        if state.marked.base[row, col] != int(2):
            col = -1
        count += 1
        state.path[count, 0] = state.path[count - 1, 0]
        state.path[count, 1] = col

    # Convert paths
    for i in range(count + 1):
        if state.marked.base[state.path[i, 0], state.path[i, 1]] == int(1):
            state.marked[state.path[i, 0], state.path[i, 1]] = 0
        else:
            state.marked[state.path[i, 0], state.path[i, 1]] = 1

    state._clear_covers()
    # Erase all prime markings
    for i in range(state.n):
        for j in range(state.m):
            if state.marked.base[i, j] == int(2):
                state.marked[i, j] = 0

    return _step3

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _step6(_HungarianState state):
    """
    Add the value found in Step 4 to every element of each covered row,
    and subtract it from every element of each uncovered column.
    Return to Step 4 without altering any stars, primes, or covered lines.
    """
    cdef double minvals = np.max(state.C.base)
    cdef np.ndarray[double, ndim=1] minval = np.asarray([np.max(state.C.base) for _ in range(state.m)])
    tmp_C = np.asarray(state.C)

    # the smallest uncovered value in the matrix
    if np.any(state.row_uncovered) and np.any(state.col_uncovered):
        for i in range(state.n):
            if state.row_uncovered[i] == 1:
                for j in range(state.m):
                    minval[j] = min(minval[j], tmp_C[i, j])

        for i in range(state.m):
            if state.col_uncovered[i] == 1:
                minvals = min(minvals, minval[i])

        for i in range(state.n):
            if state.row_uncovered[i] == 0:
                tmp_C[i, :] += minvals

            for j in range(state.m):
                if state.col_uncovered[j] == 1:
                    tmp_C[i, j] -= minvals

        state.C = tmp_C
    return _step4
