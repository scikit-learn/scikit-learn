"""Compressed Sparse Row matrix format"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['csr_matrix', 'isspmatrix_csr']


import numpy as np
from scipy._lib.six import xrange

from .base import spmatrix

from ._sparsetools import csr_tocsc, csr_tobsr, csr_count_blocks, \
        get_csr_submatrix, csr_sample_values
from .sputils import (upcast, isintlike, IndexMixin, issequence,
                      get_index_dtype, ismatrix)

from .compressed import _cs_matrix


class csr_matrix(_cs_matrix, IndexMixin):
    """
    Compressed Sparse Row matrix

    This can be instantiated in several ways:
        csr_matrix(D)
            with a dense matrix or rank-2 ndarray D

        csr_matrix(S)
            with another sparse matrix S (equivalent to S.tocsr())

        csr_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        csr_matrix((data, (row_ind, col_ind)), [shape=(M, N)])
            where ``data``, ``row_ind`` and ``col_ind`` satisfy the
            relationship ``a[row_ind[k], col_ind[k]] = data[k]``.

        csr_matrix((data, indices, indptr), [shape=(M, N)])
            is the standard CSR representation where the column indices for
            row i are stored in ``indices[indptr[i]:indptr[i+1]]`` and their
            corresponding values are stored in ``data[indptr[i]:indptr[i+1]]``.
            If the shape parameter is not supplied, the matrix dimensions
            are inferred from the index arrays.

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
        CSR format data array of the matrix
    indices
        CSR format index array of the matrix
    indptr
        CSR format index pointer array of the matrix
    has_sorted_indices
        Whether indices are sorted

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Advantages of the CSR format
      - efficient arithmetic operations CSR + CSR, CSR * CSR, etc.
      - efficient row slicing
      - fast matrix vector products

    Disadvantages of the CSR format
      - slow column slicing operations (consider CSC)
      - changes to the sparsity structure are expensive (consider LIL or DOK)

    Examples
    --------

    >>> import numpy as np
    >>> from scipy.sparse import csr_matrix
    >>> csr_matrix((3, 4), dtype=np.int8).toarray()
    array([[0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]], dtype=int8)

    >>> row = np.array([0, 0, 1, 2, 2, 2])
    >>> col = np.array([0, 2, 2, 0, 1, 2])
    >>> data = np.array([1, 2, 3, 4, 5, 6])
    >>> csr_matrix((data, (row, col)), shape=(3, 3)).toarray()
    array([[1, 0, 2],
           [0, 0, 3],
           [4, 5, 6]])

    >>> indptr = np.array([0, 2, 3, 6])
    >>> indices = np.array([0, 2, 2, 0, 1, 2])
    >>> data = np.array([1, 2, 3, 4, 5, 6])
    >>> csr_matrix((data, indices, indptr), shape=(3, 3)).toarray()
    array([[1, 0, 2],
           [0, 0, 3],
           [4, 5, 6]])

    As an example of how to construct a CSR matrix incrementally,
    the following snippet builds a term-document matrix from texts:

    >>> docs = [["hello", "world", "hello"], ["goodbye", "cruel", "world"]]
    >>> indptr = [0]
    >>> indices = []
    >>> data = []
    >>> vocabulary = {}
    >>> for d in docs:
    ...     for term in d:
    ...         index = vocabulary.setdefault(term, len(vocabulary))
    ...         indices.append(index)
    ...         data.append(1)
    ...     indptr.append(len(indices))
    ...
    >>> csr_matrix((data, indices, indptr), dtype=int).toarray()
    array([[2, 1, 0, 0],
           [0, 1, 1, 1]])

    """
    format = 'csr'

    def transpose(self, axes=None, copy=False):
        if axes is not None:
            raise ValueError(("Sparse matrices do not support "
                              "an 'axes' parameter because swapping "
                              "dimensions is the only logical permutation."))

        M, N = self.shape

        from .csc import csc_matrix
        return csc_matrix((self.data, self.indices,
                           self.indptr), shape=(N, M), copy=copy)

    transpose.__doc__ = spmatrix.transpose.__doc__

    def tolil(self, copy=False):
        from .lil import lil_matrix
        lil = lil_matrix(self.shape,dtype=self.dtype)

        self.sum_duplicates()
        ptr,ind,dat = self.indptr,self.indices,self.data
        rows, data = lil.rows, lil.data

        for n in xrange(self.shape[0]):
            start = ptr[n]
            end = ptr[n+1]
            rows[n] = ind[start:end].tolist()
            data[n] = dat[start:end].tolist()

        return lil

    tolil.__doc__ = spmatrix.tolil.__doc__

    def tocsr(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    tocsr.__doc__ = spmatrix.tocsr.__doc__

    def tocsc(self, copy=False):
        idx_dtype = get_index_dtype((self.indptr, self.indices),
                                    maxval=max(self.nnz, self.shape[0]))
        indptr = np.empty(self.shape[1] + 1, dtype=idx_dtype)
        indices = np.empty(self.nnz, dtype=idx_dtype)
        data = np.empty(self.nnz, dtype=upcast(self.dtype))

        csr_tocsc(self.shape[0], self.shape[1],
                  self.indptr.astype(idx_dtype),
                  self.indices.astype(idx_dtype),
                  self.data,
                  indptr,
                  indices,
                  data)

        from .csc import csc_matrix
        A = csc_matrix((data, indices, indptr), shape=self.shape)
        A.has_sorted_indices = True
        return A

    tocsc.__doc__ = spmatrix.tocsc.__doc__

    def tobsr(self, blocksize=None, copy=True):
        from .bsr import bsr_matrix

        if blocksize is None:
            from .spfuncs import estimate_blocksize
            return self.tobsr(blocksize=estimate_blocksize(self))

        elif blocksize == (1,1):
            arg1 = (self.data.reshape(-1,1,1),self.indices,self.indptr)
            return bsr_matrix(arg1, shape=self.shape, copy=copy)

        else:
            R,C = blocksize
            M,N = self.shape

            if R < 1 or C < 1 or M % R != 0 or N % C != 0:
                raise ValueError('invalid blocksize %s' % blocksize)

            blks = csr_count_blocks(M,N,R,C,self.indptr,self.indices)

            idx_dtype = get_index_dtype((self.indptr, self.indices),
                                        maxval=max(N//C, blks))
            indptr = np.empty(M//R+1, dtype=idx_dtype)
            indices = np.empty(blks, dtype=idx_dtype)
            data = np.zeros((blks,R,C), dtype=self.dtype)

            csr_tobsr(M, N, R, C,
                      self.indptr.astype(idx_dtype),
                      self.indices.astype(idx_dtype),
                      self.data,
                      indptr, indices, data.ravel())

            return bsr_matrix((data,indices,indptr), shape=self.shape)

    tobsr.__doc__ = spmatrix.tobsr.__doc__

    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self, x):
        """swap the members of x if this is a column-oriented matrix
        """
        return x

    def __getitem__(self, key):
        def asindices(x):
            try:
                x = np.asarray(x)

                # Check index contents to avoid creating 64bit arrays needlessly
                idx_dtype = get_index_dtype((x,), check_contents=True)
                if idx_dtype != x.dtype:
                    x = x.astype(idx_dtype)
            except:
                raise IndexError('invalid index')
            else:
                return x

        def check_bounds(indices, N):
            if indices.size == 0:
                return (0, 0)

            max_indx = indices.max()
            if max_indx >= N:
                raise IndexError('index (%d) out of range' % max_indx)

            min_indx = indices.min()
            if min_indx < -N:
                raise IndexError('index (%d) out of range' % (N + min_indx))

            return min_indx, max_indx

        def extractor(indices,N):
            """Return a sparse matrix P so that P*self implements
            slicing of the form self[[1,2,3],:]
            """
            indices = asindices(indices).copy()

            min_indx, max_indx = check_bounds(indices, N)

            if min_indx < 0:
                indices[indices < 0] += N

            indptr = np.arange(len(indices)+1, dtype=indices.dtype)
            data = np.ones(len(indices), dtype=self.dtype)
            shape = (len(indices),N)

            return csr_matrix((data,indices,indptr), shape=shape,
                              dtype=self.dtype, copy=False)

        row, col = self._unpack_index(key)

        # First attempt to use original row optimized methods
        # [1, ?]
        if isintlike(row):
            # [i, j]
            if isintlike(col):
                return self._get_single_element(row, col)
            # [i, 1:2]
            elif isinstance(col, slice):
                return self._get_row_slice(row, col)
            # [i, [1, 2]]
            elif issequence(col):
                P = extractor(col,self.shape[1]).T
                return self[row, :] * P
        elif isinstance(row, slice):
            # [1:2,??]
            if ((isintlike(col) and row.step in (1, None)) or
                    (isinstance(col, slice) and
                     col.step in (1, None) and
                     row.step in (1, None))):
                # col is int or slice with step 1, row is slice with step 1.
                return self._get_submatrix(row, col)
            elif issequence(col):
                # row is slice, col is sequence.
                P = extractor(col,self.shape[1]).T        # [1:2,[1,2]]
                sliced = self
                if row != slice(None, None, None):
                    sliced = sliced[row,:]
                return sliced * P

        elif issequence(row):
            # [[1,2],??]
            if isintlike(col) or isinstance(col,slice):
                P = extractor(row, self.shape[0])     # [[1,2],j] or [[1,2],1:2]
                extracted = P * self
                if col == slice(None, None, None):
                    return extracted
                else:
                    return extracted[:,col]

        elif ismatrix(row) and issequence(col):
            if len(row[0]) == 1 and isintlike(row[0][0]):
                # [[[1],[2]], [1,2]], outer indexing
                row = asindices(row)
                P_row = extractor(row[:,0], self.shape[0])
                P_col = extractor(col, self.shape[1]).T
                return P_row * self * P_col

        if not (issequence(col) and issequence(row)):
            # Sample elementwise
            row, col = self._index_to_arrays(row, col)

        row = asindices(row)
        col = asindices(col)
        if row.shape != col.shape:
            raise IndexError('number of row and column indices differ')
        assert row.ndim <= 2

        num_samples = np.size(row)
        if num_samples == 0:
            return csr_matrix(np.atleast_2d(row).shape, dtype=self.dtype)
        check_bounds(row, self.shape[0])
        check_bounds(col, self.shape[1])

        val = np.empty(num_samples, dtype=self.dtype)
        csr_sample_values(self.shape[0], self.shape[1],
                          self.indptr, self.indices, self.data,
                          num_samples, row.ravel(), col.ravel(), val)
        if row.ndim == 1:
            # row and col are 1d
            return np.asmatrix(val)
        return self.__class__(val.reshape(row.shape))

    def __iter__(self):
        indptr = np.zeros(2, dtype=self.indptr.dtype)
        shape = (1, self.shape[1])
        i0 = 0
        for i1 in self.indptr[1:]:
            indptr[1] = i1 - i0
            indices = self.indices[i0:i1]
            data = self.data[i0:i1]
            yield csr_matrix((data, indices, indptr), shape=shape, copy=True)
            i0 = i1

    def getrow(self, i):
        """Returns a copy of row i of the matrix, as a (1 x n)
        CSR matrix (row vector).
        """
        M, N = self.shape
        i = int(i)
        if i < 0:
            i += M
        if i < 0 or i >= M:
            raise IndexError('index (%d) out of range' % i)
        idx = slice(*self.indptr[i:i+2])
        data = self.data[idx].copy()
        indices = self.indices[idx].copy()
        indptr = np.array([0, len(indices)], dtype=self.indptr.dtype)
        return csr_matrix((data, indices, indptr), shape=(1, N),
                          dtype=self.dtype, copy=False)

    def getcol(self, i):
        """Returns a copy of column i of the matrix, as a (m x 1)
        CSR matrix (column vector).
        """
        return self._get_submatrix(slice(None), i)

    def _get_row_slice(self, i, cslice):
        """Returns a copy of row self[i, cslice]
        """
        M, N = self.shape

        if i < 0:
            i += M

        if i < 0 or i >= M:
            raise IndexError('index (%d) out of range' % i)

        start, stop, stride = cslice.indices(N)

        if stride == 1:
            # for stride == 1, get_csr_submatrix is faster
            row_indptr, row_indices, row_data = get_csr_submatrix(
                M, N, self.indptr, self.indices, self.data, i, i + 1,
                start, stop)
        else:
            # other strides need new code
            row_indices = self.indices[self.indptr[i]:self.indptr[i + 1]]
            row_data = self.data[self.indptr[i]:self.indptr[i + 1]]

            if stride > 0:
                ind = (row_indices >= start) & (row_indices < stop)
            else:
                ind = (row_indices <= start) & (row_indices > stop)

            if abs(stride) > 1:
                ind &= (row_indices - start) % stride == 0

            row_indices = (row_indices[ind] - start) // stride
            row_data = row_data[ind]
            row_indptr = np.array([0, len(row_indices)])

            if stride < 0:
                row_data = row_data[::-1]
                row_indices = abs(row_indices[::-1])

        shape = (1, int(np.ceil(float(stop - start) / stride)))
        return csr_matrix((row_data, row_indices, row_indptr), shape=shape,
                          dtype=self.dtype, copy=False)

    def _get_submatrix(self, row_slice, col_slice):
        """Return a submatrix of this matrix (new matrix is created)."""

        def process_slice(sl, num):
            if isinstance(sl, slice):
                i0, i1, stride = sl.indices(num)
                if stride != 1:
                    raise ValueError('slicing with step != 1 not supported')
            elif isintlike(sl):
                if sl < 0:
                    sl += num
                i0, i1 = sl, sl + 1
            else:
                raise TypeError('expected slice or scalar')

            if not (0 <= i0 <= num) or not (0 <= i1 <= num) or not (i0 <= i1):
                raise IndexError(
                      "index out of bounds: 0 <= %d <= %d, 0 <= %d <= %d,"
                      " %d <= %d" % (i0, num, i1, num, i0, i1))
            return i0, i1

        M,N = self.shape
        i0, i1 = process_slice(row_slice, M)
        j0, j1 = process_slice(col_slice, N)

        indptr, indices, data = get_csr_submatrix(
            M, N, self.indptr, self.indices, self.data, i0, i1, j0, j1)

        shape = (i1 - i0, j1 - j0)
        return self.__class__((data, indices, indptr), shape=shape,
                              dtype=self.dtype, copy=False)


def isspmatrix_csr(x):
    """Is x of csr_matrix type?

    Parameters
    ----------
    x
        object to check for being a csr matrix

    Returns
    -------
    bool
        True if x is a csr matrix, False otherwise

    Examples
    --------
    >>> from scipy.sparse import csr_matrix, isspmatrix_csr
    >>> isspmatrix_csr(csr_matrix([[5]]))
    True

    >>> from scipy.sparse import csc_matrix, csr_matrix, isspmatrix_csc
    >>> isspmatrix_csr(csc_matrix([[5]]))
    False
    """
    return isinstance(x, csr_matrix)
