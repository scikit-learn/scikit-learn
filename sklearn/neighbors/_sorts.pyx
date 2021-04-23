# distutils : language = c++

# Some internals rely on some kinds of sorts, like KNeighborsMixin for
# partitioning neighbors.
#
# The C++ std library exposes a lot of efficient sorting algorithms,
# like nth_element, which is an efficient partial sort for KNeighborsMixin's
# case.
#
# To use std::algorithm, a few fixture can be defined using Cython, mainly:
# - Cython functions used in Cython implementations. Those call
# - C++ functions that wraps function of std::algorithm and (sometimes) use
# - an Comparator to state how to sort
#
# We could directly call the C++ interfaces (*_inner) in cython.
# This works but we prefer to have the cython interface declared in
# headers files (*.pxd) and used in Cython code, especially to hide
# C++ templating from Cython callers.
#

import numpy as np
cimport numpy as np
from cython.parallel import prange, parallel

from cython cimport floating, integral

cdef extern from *:
    """
    #include <algorithm>

    template<class D, class I>
    class IndexComparator {
    private:
        const D *data;
    public:
        IndexComparator(const D *data):
            data(data) {}

        bool operator()(const I &a, const I &b) const {
            return data[a] == data[b] ? a < b : data[a] < data[b];
        }
    };

    template<class D, class I>
    void intro_select_inner(
        const D *data,
        I *indices,
        const I &pivot,
        const I &n_points) {
        IndexComparator<D, I> index_comparator(data);
        std::nth_element(
            indices,
            indices + pivot,
            indices + n_points,
            index_comparator);
    }
    """
    void intro_select_inner[D, I](
                D *data,
                I *indices,
                I pivot,
                I n_points) nogil except +


cdef integral intro_select(
        floating *data,
        integral *indices,
        integral pivot,
        integral n_points) nogil:
    """Partition indices based on their associated data.

    It is essentially a partial in-place quicksort around a
    set pivot, that is upon return, the values in indices will
    be rearranged such that (assuming numpy-style indexing):

    data[indices[:pivot]] <= data[indices[pivot]] <= data[indices[pivot:]]

    Parameters
    ----------
    data : floating pointer
        Pointer to a 1D array of length n_points containing floating data
    indices : int pointer
        Pointer to a 1D array of length n_points. This will be modified
        in-place.
    pivot : int
        the index within indices around which to split the points.
    n_points : int
        the length of data and indices.
    Returns
    -------
    status : int
        integer exit status.  On return, the contents of indices are
        modified as noted above.
    """
    intro_select_inner(
        data,
        indices,
        pivot,
        n_points,
    )
    return 0

cpdef np.ndarray[integral, ndim=2, mode='c'] argpartition(
        np.ndarray[floating, ndim=2, mode='c'] data,
        integral pivot):
    """
    Return an array of indices such that selection of these indices on
    the original data would return a row-wise partitioned version of
    the array with respect to the value of the pivot.
    
    This is equivalent to using:
        
        np.argpartition(data, pivot, axis=1)
    
    Notes
    -----
    Like np.argpartition, this also makes use of intro_select
    but via the implementation in std::algorithm (nth_element).
    
    Hence resulting indices might be different.
    
    Parameters
    ----------
    data : array of floating
        A 2D array containing floating data.
    pivot : int
        The index around which to sort row indices.

    Returns
    -------
    indices : array of integer
         A 2D array of row indices partitioned with respect to the given pivot.
    """
    cdef:
        integral i_row
        integral n_rows = data.shape[0]
        integral n_cols = data.shape[1]
        floating *data_ptr = &data[0, 0]

    cdef np.ndarray[integral, ndim=2, mode='c'] indices = np.tile(
        np.arange(n_cols, dtype=int), reps=(n_rows, 1)
    )
    cdef integral * indices_ptr = &indices[0, 0]

    # Sorting on rows in parallel, indices is modified inplace
    for i_row in prange(n_rows, schedule='static', nogil=True):
        intro_select(
            data_ptr + i_row * n_cols,
            indices_ptr + i_row * n_cols,
            pivot,
            n_cols,
        )
    return indices


