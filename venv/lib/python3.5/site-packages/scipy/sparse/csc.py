"""Compressed Sparse Column matrix format"""
from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['csc_matrix', 'isspmatrix_csc']


import numpy as np

from .base import spmatrix
from ._sparsetools import csc_tocsr
from . import _sparsetools
from .sputils import upcast, isintlike, IndexMixin, get_index_dtype

from .compressed import _cs_matrix


class csc_matrix(_cs_matrix, IndexMixin):
    """
    Compressed Sparse Column matrix

    This can be instantiated in several ways:

        csc_matrix(D)
            with a dense matrix or rank-2 ndarray D

        csc_matrix(S)
            with another sparse matrix S (equivalent to S.tocsc())

        csc_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        csc_matrix((data, (row_ind, col_ind)), [shape=(M, N)])
            where ``data``, ``row_ind`` and ``col_ind`` satisfy the
            relationship ``a[row_ind[k], col_ind[k]] = data[k]``.

        csc_matrix((data, indices, indptr), [shape=(M, N)])
            is the standard CSC representation where the row indices for
            column i are stored in ``indices[indptr[i]:indptr[i+1]]``
            and their corresponding values are stored in
            ``data[indptr[i]:indptr[i+1]]``.  If the shape parameter is
            not supplied, the matrix dimensions are inferred from
            the index arrays.

    Attributes
    ----------
    dtype : dtype
        Data type of the matrix
    shape : 2-tuple
        Shape of the matrix
    ndim : int
        Number of dimensions (this is always 2)
    nnz
        Number of nonzero elements
    data
        Data array of the matrix
    indices
        CSC format index array
    indptr
        CSC format index pointer array
    has_sorted_indices
        Whether indices are sorted

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Advantages of the CSC format
        - efficient arithmetic operations CSC + CSC, CSC * CSC, etc.
        - efficient column slicing
        - fast matrix vector products (CSR, BSR may be faster)

    Disadvantages of the CSC format
      - slow row slicing operations (consider CSR)
      - changes to the sparsity structure are expensive (consider LIL or DOK)


    Examples
    --------

    >>> import numpy as np
    >>> from scipy.sparse import csc_matrix
    >>> csc_matrix((3, 4), dtype=np.int8).toarray()
    array([[0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]], dtype=int8)

    >>> row = np.array([0, 2, 2, 0, 1, 2])
    >>> col = np.array([0, 0, 1, 2, 2, 2])
    >>> data = np.array([1, 2, 3, 4, 5, 6])
    >>> csc_matrix((data, (row, col)), shape=(3, 3)).toarray()
    array([[1, 0, 4],
           [0, 0, 5],
           [2, 3, 6]])

    >>> indptr = np.array([0, 2, 3, 6])
    >>> indices = np.array([0, 2, 2, 0, 1, 2])
    >>> data = np.array([1, 2, 3, 4, 5, 6])
    >>> csc_matrix((data, indices, indptr), shape=(3, 3)).toarray()
    array([[1, 0, 4],
           [0, 0, 5],
           [2, 3, 6]])

    """
    format = 'csc'

    def transpose(self, axes=None, copy=False):
        if axes is not None:
            raise ValueError(("Sparse matrices do not support "
                              "an 'axes' parameter because swapping "
                              "dimensions is the only logical permutation."))

        M, N = self.shape

        from .csr import csr_matrix
        return csr_matrix((self.data, self.indices,
                           self.indptr), (N, M), copy=copy)

    transpose.__doc__ = spmatrix.transpose.__doc__

    def __iter__(self):
        for r in self.tocsr():
            yield r

    def tocsc(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    tocsc.__doc__ = spmatrix.tocsc.__doc__

    def tocsr(self, copy=False):
        M,N = self.shape
        idx_dtype = get_index_dtype((self.indptr, self.indices),
                                    maxval=max(self.nnz, N))
        indptr = np.empty(M + 1, dtype=idx_dtype)
        indices = np.empty(self.nnz, dtype=idx_dtype)
        data = np.empty(self.nnz, dtype=upcast(self.dtype))

        csc_tocsr(M, N,
                  self.indptr.astype(idx_dtype),
                  self.indices.astype(idx_dtype),
                  self.data,
                  indptr,
                  indices,
                  data)

        from .csr import csr_matrix
        A = csr_matrix((data, indices, indptr), shape=self.shape, copy=False)
        A.has_sorted_indices = True
        return A

    tocsr.__doc__ = spmatrix.tocsr.__doc__

    def __getitem__(self, key):
        # Use CSR to implement fancy indexing.

        row, col = self._unpack_index(key)
        # Things that return submatrices. row or col is a int or slice.
        if (isinstance(row, slice) or isinstance(col, slice) or
                isintlike(row) or isintlike(col)):
            return self.T[col, row].T
        # Things that return a sequence of values.
        else:
            return self.T[col, row]

    def nonzero(self):
        # CSC can't use _cs_matrix's .nonzero method because it
        # returns the indices sorted for self transposed.

        # Get row and col indices, from _cs_matrix.tocoo
        major_dim, minor_dim = self._swap(self.shape)
        minor_indices = self.indices
        major_indices = np.empty(len(minor_indices), dtype=self.indices.dtype)
        _sparsetools.expandptr(major_dim, self.indptr, major_indices)
        row, col = self._swap((major_indices, minor_indices))

        # Remove explicit zeros
        nz_mask = self.data != 0
        row = row[nz_mask]
        col = col[nz_mask]

        # Sort them to be in C-style order
        ind = np.argsort(row, kind='mergesort')
        row = row[ind]
        col = col[ind]

        return row, col

    nonzero.__doc__ = _cs_matrix.nonzero.__doc__

    def getrow(self, i):
        """Returns a copy of row i of the matrix, as a (1 x n)
        CSR matrix (row vector).
        """
        # we convert to CSR to maintain compatibility with old impl.
        # in spmatrix.getrow()
        return self._get_submatrix(i, slice(None)).tocsr()

    def getcol(self, i):
        """Returns a copy of column i of the matrix, as a (m x 1)
        CSC matrix (column vector).
        """
        M, N = self.shape
        i = int(i)
        if i < 0:
            i += N
        if i < 0 or i >= N:
            raise IndexError('index (%d) out of range' % i)
        idx = slice(*self.indptr[i:i+2])
        data = self.data[idx].copy()
        indices = self.indices[idx].copy()
        indptr = np.array([0, len(indices)], dtype=self.indptr.dtype)
        return csc_matrix((data, indices, indptr), shape=(M, 1),
                          dtype=self.dtype, copy=False)

    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self, x):
        """swap the members of x if this is a column-oriented matrix
        """
        return x[1], x[0]


def isspmatrix_csc(x):
    """Is x of csc_matrix type?

    Parameters
    ----------
    x
        object to check for being a csc matrix

    Returns
    -------
    bool
        True if x is a csc matrix, False otherwise

    Examples
    --------
    >>> from scipy.sparse import csc_matrix, isspmatrix_csc
    >>> isspmatrix_csc(csc_matrix([[5]]))
    True

    >>> from scipy.sparse import csc_matrix, csr_matrix, isspmatrix_csc
    >>> isspmatrix_csc(csr_matrix([[5]]))
    False
    """
    return isinstance(x, csc_matrix)
