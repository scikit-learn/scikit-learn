"""Sparse DIAgonal format"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['dia_matrix', 'isspmatrix_dia']

import numpy as np

from .base import isspmatrix, _formats, spmatrix
from .data import _data_matrix
from .sputils import (isshape, upcast_char, getdtype, get_index_dtype,
                      get_sum_dtype, validateaxis, check_shape)
from ._sparsetools import dia_matvec


class dia_matrix(_data_matrix):
    """Sparse matrix with DIAgonal storage

    This can be instantiated in several ways:
        dia_matrix(D)
            with a dense matrix

        dia_matrix(S)
            with another sparse matrix S (equivalent to S.todia())

        dia_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N),
            dtype is optional, defaulting to dtype='d'.

        dia_matrix((data, offsets), shape=(M, N))
            where the ``data[k,:]`` stores the diagonal entries for
            diagonal ``offsets[k]`` (See example below)

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
        DIA format data array of the matrix
    offsets
        DIA format offset array of the matrix

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Examples
    --------

    >>> import numpy as np
    >>> from scipy.sparse import dia_matrix
    >>> dia_matrix((3, 4), dtype=np.int8).toarray()
    array([[0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]], dtype=int8)

    >>> data = np.array([[1, 2, 3, 4]]).repeat(3, axis=0)
    >>> offsets = np.array([0, -1, 2])
    >>> dia_matrix((data, offsets), shape=(4, 4)).toarray()
    array([[1, 0, 3, 0],
           [1, 2, 0, 4],
           [0, 2, 3, 0],
           [0, 0, 3, 4]])

    """
    format = 'dia'

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        _data_matrix.__init__(self)

        if isspmatrix_dia(arg1):
            if copy:
                arg1 = arg1.copy()
            self.data = arg1.data
            self.offsets = arg1.offsets
            self._shape = check_shape(arg1.shape)
        elif isspmatrix(arg1):
            if isspmatrix_dia(arg1) and copy:
                A = arg1.copy()
            else:
                A = arg1.todia()
            self.data = A.data
            self.offsets = A.offsets
            self._shape = check_shape(A.shape)
        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                # create empty matrix
                self._shape = check_shape(arg1)
                self.data = np.zeros((0,0), getdtype(dtype, default=float))
                idx_dtype = get_index_dtype(maxval=max(self.shape))
                self.offsets = np.zeros((0), dtype=idx_dtype)
            else:
                try:
                    # Try interpreting it as (data, offsets)
                    data, offsets = arg1
                except:
                    raise ValueError('unrecognized form for dia_matrix constructor')
                else:
                    if shape is None:
                        raise ValueError('expected a shape argument')
                    self.data = np.atleast_2d(np.array(arg1[0], dtype=dtype, copy=copy))
                    self.offsets = np.atleast_1d(np.array(arg1[1],
                                                          dtype=get_index_dtype(maxval=max(shape)),
                                                          copy=copy))
                    self._shape = check_shape(shape)
        else:
            #must be dense, convert to COO first, then to DIA
            try:
                arg1 = np.asarray(arg1)
            except:
                raise ValueError("unrecognized form for"
                        " %s_matrix constructor" % self.format)
            from .coo import coo_matrix
            A = coo_matrix(arg1, dtype=dtype, shape=shape).todia()
            self.data = A.data
            self.offsets = A.offsets
            self._shape = check_shape(A.shape)

        if dtype is not None:
            self.data = self.data.astype(dtype)

        #check format
        if self.offsets.ndim != 1:
            raise ValueError('offsets array must have rank 1')

        if self.data.ndim != 2:
            raise ValueError('data array must have rank 2')

        if self.data.shape[0] != len(self.offsets):
            raise ValueError('number of diagonals (%d) '
                    'does not match the number of offsets (%d)'
                    % (self.data.shape[0], len(self.offsets)))

        if len(np.unique(self.offsets)) != len(self.offsets):
            raise ValueError('offset array contains duplicate values')

    def __repr__(self):
        format = _formats[self.getformat()][1]
        return "<%dx%d sparse matrix of type '%s'\n" \
               "\twith %d stored elements (%d diagonals) in %s format>" % \
               (self.shape + (self.dtype.type, self.nnz, self.data.shape[0],
                              format))

    def _data_mask(self):
        """Returns a mask of the same shape as self.data, where
        mask[i,j] is True when data[i,j] corresponds to a stored element."""
        num_rows, num_cols = self.shape
        offset_inds = np.arange(self.data.shape[1])
        row = offset_inds - self.offsets[:,None]
        mask = (row >= 0)
        mask &= (row < num_rows)
        mask &= (offset_inds < num_cols)
        return mask

    def count_nonzero(self):
        mask = self._data_mask()
        return np.count_nonzero(self.data[mask])

    def getnnz(self, axis=None):
        if axis is not None:
            raise NotImplementedError("getnnz over an axis is not implemented "
                                      "for DIA format")
        M,N = self.shape
        nnz = 0
        for k in self.offsets:
            if k > 0:
                nnz += min(M,N-k)
            else:
                nnz += min(M+k,N)
        return int(nnz)

    getnnz.__doc__ = spmatrix.getnnz.__doc__
    count_nonzero.__doc__ = spmatrix.count_nonzero.__doc__

    def sum(self, axis=None, dtype=None, out=None):
        validateaxis(axis)

        if axis is not None and axis < 0:
            axis += 2

        res_dtype = get_sum_dtype(self.dtype)
        num_rows, num_cols = self.shape
        ret = None

        if axis == 0:
            mask = self._data_mask()
            x = (self.data * mask).sum(axis=0)
            if x.shape[0] == num_cols:
                res = x
            else:
                res = np.zeros(num_cols, dtype=x.dtype)
                res[:x.shape[0]] = x
            ret = np.matrix(res, dtype=res_dtype)

        else:
            row_sums = np.zeros(num_rows, dtype=res_dtype)
            one = np.ones(num_cols, dtype=res_dtype)
            dia_matvec(num_rows, num_cols, len(self.offsets),
                       self.data.shape[1], self.offsets, self.data, one, row_sums)

            row_sums = np.matrix(row_sums)

            if axis is None:
                return row_sums.sum(dtype=dtype, out=out)

            if axis is not None:
                row_sums = row_sums.T

            ret = np.matrix(row_sums.sum(axis=axis))

        if out is not None and out.shape != ret.shape:
            raise ValueError("dimensions do not match")

        return ret.sum(axis=(), dtype=dtype, out=out)

    sum.__doc__ = spmatrix.sum.__doc__

    def _mul_vector(self, other):
        x = other

        y = np.zeros(self.shape[0], dtype=upcast_char(self.dtype.char,
                                                       x.dtype.char))

        L = self.data.shape[1]

        M,N = self.shape

        dia_matvec(M,N, len(self.offsets), L, self.offsets, self.data, x.ravel(), y.ravel())

        return y

    def _mul_multimatrix(self, other):
        return np.hstack([self._mul_vector(col).reshape(-1,1) for col in other.T])

    def _setdiag(self, values, k=0):
        M, N = self.shape

        if values.ndim == 0:
            # broadcast
            values_n = np.inf
        else:
            values_n = len(values)

        if k < 0:
            n = min(M + k, N, values_n)
            min_index = 0
            max_index = n
        else:
            n = min(M, N - k, values_n)
            min_index = k
            max_index = k + n

        if values.ndim != 0:
            # allow also longer sequences
            values = values[:n]

        if k in self.offsets:
            self.data[self.offsets == k, min_index:max_index] = values
        else:
            self.offsets = np.append(self.offsets, self.offsets.dtype.type(k))
            m = max(max_index, self.data.shape[1])
            data = np.zeros((self.data.shape[0]+1, m), dtype=self.data.dtype)
            data[:-1,:self.data.shape[1]] = self.data
            data[-1, min_index:max_index] = values
            self.data = data

    def todia(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    todia.__doc__ = spmatrix.todia.__doc__

    def transpose(self, axes=None, copy=False):
        if axes is not None:
            raise ValueError(("Sparse matrices do not support "
                              "an 'axes' parameter because swapping "
                              "dimensions is the only logical permutation."))

        num_rows, num_cols = self.shape
        max_dim = max(self.shape)

        # flip diagonal offsets
        offsets = -self.offsets

        # re-align the data matrix
        r = np.arange(len(offsets), dtype=np.intc)[:, None]
        c = np.arange(num_rows, dtype=np.intc) - (offsets % max_dim)[:, None]
        pad_amount = max(0, max_dim-self.data.shape[1])
        data = np.hstack((self.data, np.zeros((self.data.shape[0], pad_amount),
                                              dtype=self.data.dtype)))
        data = data[r, c]
        return dia_matrix((data, offsets), shape=(
            num_cols, num_rows), copy=copy)

    transpose.__doc__ = spmatrix.transpose.__doc__

    def diagonal(self, k=0):
        rows, cols = self.shape
        if k <= -rows or k >= cols:
            raise ValueError("k exceeds matrix dimensions")
        idx, = np.where(self.offsets == k)
        first_col, last_col = max(0, k), min(rows + k, cols)
        if idx.size == 0:
            return np.zeros(last_col - first_col, dtype=self.data.dtype)
        return self.data[idx[0], first_col:last_col]

    diagonal.__doc__ = spmatrix.diagonal.__doc__

    def tocsc(self, copy=False):
        from .csc import csc_matrix
        if self.nnz == 0:
            return csc_matrix(self.shape, dtype=self.dtype)

        num_rows, num_cols = self.shape
        num_offsets, offset_len = self.data.shape
        offset_inds = np.arange(offset_len)

        row = offset_inds - self.offsets[:,None]
        mask = (row >= 0)
        mask &= (row < num_rows)
        mask &= (offset_inds < num_cols)
        mask &= (self.data != 0)

        idx_dtype = get_index_dtype(maxval=max(self.shape))
        indptr = np.zeros(num_cols + 1, dtype=idx_dtype)
        indptr[1:offset_len+1] = np.cumsum(mask.sum(axis=0))
        indptr[offset_len+1:] = indptr[offset_len]
        indices = row.T[mask.T].astype(idx_dtype, copy=False)
        data = self.data.T[mask.T]
        return csc_matrix((data, indices, indptr), shape=self.shape,
                          dtype=self.dtype)

    tocsc.__doc__ = spmatrix.tocsc.__doc__

    def tocoo(self, copy=False):
        num_rows, num_cols = self.shape
        num_offsets, offset_len = self.data.shape
        offset_inds = np.arange(offset_len)

        row = offset_inds - self.offsets[:,None]
        mask = (row >= 0)
        mask &= (row < num_rows)
        mask &= (offset_inds < num_cols)
        mask &= (self.data != 0)
        row = row[mask]
        col = np.tile(offset_inds, num_offsets)[mask.ravel()]
        data = self.data[mask]

        from .coo import coo_matrix
        A = coo_matrix((data,(row,col)), shape=self.shape, dtype=self.dtype)
        A.has_canonical_format = True
        return A

    tocoo.__doc__ = spmatrix.tocoo.__doc__

    # needed by _data_matrix
    def _with_data(self, data, copy=True):
        """Returns a matrix with the same sparsity structure as self,
        but with different data.  By default the structure arrays are copied.
        """
        if copy:
            return dia_matrix((data, self.offsets.copy()), shape=self.shape)
        else:
            return dia_matrix((data,self.offsets), shape=self.shape)

    def resize(self, *shape):
        shape = check_shape(shape)
        M, N = shape
        # we do not need to handle the case of expanding N
        self.data = self.data[:, :N]

        if (M > self.shape[0] and
                np.any(self.offsets + self.shape[0] < self.data.shape[1])):
            # explicitly clear values that were previously hidden
            mask = (self.offsets[:, None] + self.shape[0] <=
                    np.arange(self.data.shape[1]))
            self.data[mask] = 0

        self._shape = shape

    resize.__doc__ = spmatrix.resize.__doc__


def isspmatrix_dia(x):
    """Is x of dia_matrix type?

    Parameters
    ----------
    x
        object to check for being a dia matrix

    Returns
    -------
    bool
        True if x is a dia matrix, False otherwise

    Examples
    --------
    >>> from scipy.sparse import dia_matrix, isspmatrix_dia
    >>> isspmatrix_dia(dia_matrix([[5]]))
    True

    >>> from scipy.sparse import dia_matrix, csr_matrix, isspmatrix_dia
    >>> isspmatrix_dia(csr_matrix([[5]]))
    False
    """
    return isinstance(x, dia_matrix)
