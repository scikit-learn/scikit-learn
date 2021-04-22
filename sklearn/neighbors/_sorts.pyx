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
import numpy as np
cimport numpy as np
from cython.parallel import prange, parallel

from ._typedefs import ITYPE

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


cdef int intro_select(
        DTYPE_t *data,
        ITYPE_t *indices,
        ITYPE_t pivot,
        ITYPE_t n_points) nogil except -1:
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
        n_points)
    return 0

cpdef np.ndarray[ITYPE_t, ndim=2, mode='c'] argpartition(
        np.ndarray[DTYPE_t, ndim=2, mode='c'] data,
        ITYPE_t pivot):
    cdef:
        ITYPE_t i_row
        ITYPE_t n_rows = data.shape[0]
        ITYPE_t n_cols = data.shape[1]
        DTYPE_t *data_ptr = &data[0, 0]

    # TODO: I haven't found a nicer way to create indices yet.
    cdef np.ndarray[ITYPE_t, ndim = 2, mode = 'c'] indices = (
    np.ones((n_rows, 1), dtype=ITYPE) @ np.arange(n_cols, dtype=ITYPE)[None, :]
    )
    cdef ITYPE_t * indices_ptr = &indices[0, 0]

    # Sorting on rows in parallel, indices is modified inplace
    for i_row in prange(n_rows, schedule='static', nogil=True):
        intro_select(
            data_ptr + i_row * n_cols,
            indices_ptr + i_row * n_cols,
            pivot,
            n_cols)
    return indices


