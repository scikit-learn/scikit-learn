"""Dictionary Of Keys based matrix"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['dok_matrix', 'isspmatrix_dok']

import functools
import operator
import itertools

import numpy as np

from scipy._lib.six import zip as izip, xrange, iteritems, iterkeys, itervalues

from .base import spmatrix, isspmatrix
from .sputils import (isdense, getdtype, isshape, isintlike, isscalarlike,
                      upcast, upcast_scalar, IndexMixin, get_index_dtype,
                      check_shape)

try:
    from operator import isSequenceType as _is_sequence
except ImportError:
    def _is_sequence(x):
        return (hasattr(x, '__len__') or hasattr(x, '__next__')
                or hasattr(x, 'next'))


class dok_matrix(spmatrix, IndexMixin, dict):
    """
    Dictionary Of Keys based sparse matrix.

    This is an efficient structure for constructing sparse
    matrices incrementally.

    This can be instantiated in several ways:
        dok_matrix(D)
            with a dense matrix, D

        dok_matrix(S)
            with a sparse matrix, S

        dok_matrix((M,N), [dtype])
            create the matrix with initial shape (M,N)
            dtype is optional, defaulting to dtype='d'

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

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Allows for efficient O(1) access of individual elements.
    Duplicates are not allowed.
    Can be efficiently converted to a coo_matrix once constructed.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import dok_matrix
    >>> S = dok_matrix((5, 5), dtype=np.float32)
    >>> for i in range(5):
    ...     for j in range(5):
    ...         S[i, j] = i + j    # Update element

    """
    format = 'dok'

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        dict.__init__(self)
        spmatrix.__init__(self)

        self.dtype = getdtype(dtype, default=float)
        if isinstance(arg1, tuple) and isshape(arg1):  # (M,N)
            M, N = arg1
            self._shape = check_shape((M, N))
        elif isspmatrix(arg1):  # Sparse ctor
            if isspmatrix_dok(arg1) and copy:
                arg1 = arg1.copy()
            else:
                arg1 = arg1.todok()

            if dtype is not None:
                arg1 = arg1.astype(dtype)

            dict.update(self, arg1)
            self._shape = check_shape(arg1.shape)
            self.dtype = arg1.dtype
        else:  # Dense ctor
            try:
                arg1 = np.asarray(arg1)
            except:
                raise TypeError('Invalid input format.')

            if len(arg1.shape) != 2:
                raise TypeError('Expected rank <=2 dense array or matrix.')

            from .coo import coo_matrix
            d = coo_matrix(arg1, dtype=dtype).todok()
            dict.update(self, d)
            self._shape = check_shape(arg1.shape)
            self.dtype = d.dtype

    def update(self, val):
        # Prevent direct usage of update
        raise NotImplementedError("Direct modification to dok_matrix element "
                                  "is not allowed.")

    def _update(self, data):
        """An update method for dict data defined for direct access to
        `dok_matrix` data. Main purpose is to be used for effcient conversion
        from other spmatrix classes. Has no checking if `data` is valid."""
        return dict.update(self, data)

    def set_shape(self, shape):
        new_matrix = self.reshape(shape, copy=False).asformat(self.format)
        self.__dict__ = new_matrix.__dict__
        dict.clear(self)
        dict.update(self, new_matrix)

    shape = property(fget=spmatrix.get_shape, fset=set_shape)

    def getnnz(self, axis=None):
        if axis is not None:
            raise NotImplementedError("getnnz over an axis is not implemented "
                                      "for DOK format.")
        return dict.__len__(self)

    def count_nonzero(self):
        return sum(x != 0 for x in itervalues(self))

    getnnz.__doc__ = spmatrix.getnnz.__doc__
    count_nonzero.__doc__ = spmatrix.count_nonzero.__doc__

    def __len__(self):
        return dict.__len__(self)

    def get(self, key, default=0.):
        """This overrides the dict.get method, providing type checking
        but otherwise equivalent functionality.
        """
        try:
            i, j = key
            assert isintlike(i) and isintlike(j)
        except (AssertionError, TypeError, ValueError):
            raise IndexError('Index must be a pair of integers.')
        if (i < 0 or i >= self.shape[0] or j < 0 or j >= self.shape[1]):
            raise IndexError('Index out of bounds.')
        return dict.get(self, key, default)

    def __getitem__(self, index):
        """If key=(i, j) is a pair of integers, return the corresponding
        element.  If either i or j is a slice or sequence, return a new sparse
        matrix with just these elements.
        """
        zero = self.dtype.type(0)
        i, j = self._unpack_index(index)

        i_intlike = isintlike(i)
        j_intlike = isintlike(j)

        if i_intlike and j_intlike:
            i = int(i)
            j = int(j)
            if i < 0:
                i += self.shape[0]
            if i < 0 or i >= self.shape[0]:
                raise IndexError('Index out of bounds.')
            if j < 0:
                j += self.shape[1]
            if j < 0 or j >= self.shape[1]:
                raise IndexError('Index out of bounds.')
            return dict.get(self, (i,j), zero)
        elif ((i_intlike or isinstance(i, slice)) and
              (j_intlike or isinstance(j, slice))):
            # Fast path for slicing very sparse matrices
            i_slice = slice(i, i+1) if i_intlike else i
            j_slice = slice(j, j+1) if j_intlike else j
            i_indices = i_slice.indices(self.shape[0])
            j_indices = j_slice.indices(self.shape[1])
            i_seq = xrange(*i_indices)
            j_seq = xrange(*j_indices)
            newshape = (len(i_seq), len(j_seq))
            newsize = _prod(newshape)

            if len(self) < 2*newsize and newsize != 0:
                # Switch to the fast path only when advantageous
                # (count the iterations in the loops, adjust for complexity)
                #
                # We also don't handle newsize == 0 here (if
                # i/j_intlike, it can mean index i or j was out of
                # bounds)
                return self._getitem_ranges(i_indices, j_indices, newshape)

        i, j = self._index_to_arrays(i, j)

        if i.size == 0:
            return dok_matrix(i.shape, dtype=self.dtype)

        min_i = i.min()
        if min_i < -self.shape[0] or i.max() >= self.shape[0]:
            raise IndexError('Index (%d) out of range -%d to %d.' %
                             (i.min(), self.shape[0], self.shape[0]-1))
        if min_i < 0:
            i = i.copy()
            i[i < 0] += self.shape[0]

        min_j = j.min()
        if min_j < -self.shape[1] or j.max() >= self.shape[1]:
            raise IndexError('Index (%d) out of range -%d to %d.' %
                             (j.min(), self.shape[1], self.shape[1]-1))
        if min_j < 0:
            j = j.copy()
            j[j < 0] += self.shape[1]

        newdok = dok_matrix(i.shape, dtype=self.dtype)

        for key in itertools.product(xrange(i.shape[0]), xrange(i.shape[1])):
            v = dict.get(self, (i[key], j[key]), zero)
            if v:
                dict.__setitem__(newdok, key, v)

        return newdok

    def _getitem_ranges(self, i_indices, j_indices, shape):
        # performance golf: we don't want Numpy scalars here, they are slow
        i_start, i_stop, i_stride = map(int, i_indices)
        j_start, j_stop, j_stride = map(int, j_indices)

        newdok = dok_matrix(shape, dtype=self.dtype)

        for (ii, jj) in iterkeys(self):
            # ditto for numpy scalars
            ii = int(ii)
            jj = int(jj)
            a, ra = divmod(ii - i_start, i_stride)
            if a < 0 or a >= shape[0] or ra != 0:
                continue
            b, rb = divmod(jj - j_start, j_stride)
            if b < 0 or b >= shape[1] or rb != 0:
                continue
            dict.__setitem__(newdok, (a, b),
                             dict.__getitem__(self, (ii, jj)))
        return newdok

    def __setitem__(self, index, x):
        if isinstance(index, tuple) and len(index) == 2:
            # Integer index fast path
            i, j = index
            if (isintlike(i) and isintlike(j) and 0 <= i < self.shape[0]
                    and 0 <= j < self.shape[1]):
                v = np.asarray(x, dtype=self.dtype)
                if v.ndim == 0 and v != 0:
                    dict.__setitem__(self, (int(i), int(j)), v[()])
                    return

        i, j = self._unpack_index(index)
        i, j = self._index_to_arrays(i, j)

        if isspmatrix(x):
            x = x.toarray()

        # Make x and i into the same shape
        x = np.asarray(x, dtype=self.dtype)
        x, _ = np.broadcast_arrays(x, i)

        if x.shape != i.shape:
            raise ValueError("Shape mismatch in assignment.")

        if np.size(x) == 0:
            return

        min_i = i.min()
        if min_i < -self.shape[0] or i.max() >= self.shape[0]:
            raise IndexError('Index (%d) out of range -%d to %d.' %
                             (i.min(), self.shape[0], self.shape[0]-1))
        if min_i < 0:
            i = i.copy()
            i[i < 0] += self.shape[0]

        min_j = j.min()
        if min_j < -self.shape[1] or j.max() >= self.shape[1]:
            raise IndexError('Index (%d) out of range -%d to %d.' %
                             (j.min(), self.shape[1], self.shape[1]-1))
        if min_j < 0:
            j = j.copy()
            j[j < 0] += self.shape[1]

        dict.update(self, izip(izip(i.flat, j.flat), x.flat))

        if 0 in x:
            zeroes = x == 0
            for key in izip(i[zeroes].flat, j[zeroes].flat):
                if dict.__getitem__(self, key) == 0:
                    # may have been superseded by later update
                    del self[key]

    def __add__(self, other):
        if isscalarlike(other):
            res_dtype = upcast_scalar(self.dtype, other)
            new = dok_matrix(self.shape, dtype=res_dtype)
            # Add this scalar to every element.
            M, N = self.shape
            for key in itertools.product(xrange(M), xrange(N)):
                aij = dict.get(self, (key), 0) + other
                if aij:
                    new[key] = aij
            # new.dtype.char = self.dtype.char
        elif isspmatrix_dok(other):
            if other.shape != self.shape:
                raise ValueError("Matrix dimensions are not equal.")
            # We could alternatively set the dimensions to the largest of
            # the two matrices to be summed.  Would this be a good idea?
            res_dtype = upcast(self.dtype, other.dtype)
            new = dok_matrix(self.shape, dtype=res_dtype)
            dict.update(new, self)
            with np.errstate(over='ignore'):
                dict.update(new,
                           ((k, new[k] + other[k]) for k in iterkeys(other)))
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = self.todense() + other
        else:
            return NotImplemented
        return new

    def __radd__(self, other):
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            M, N = self.shape
            for key in itertools.product(xrange(M), xrange(N)):
                aij = dict.get(self, (key), 0) + other
                if aij:
                    new[key] = aij
        elif isspmatrix_dok(other):
            if other.shape != self.shape:
                raise ValueError("Matrix dimensions are not equal.")
            new = dok_matrix(self.shape, dtype=self.dtype)
            dict.update(new, self)
            dict.update(new,
                       ((k, self[k] + other[k]) for k in iterkeys(other)))
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = other + self.todense()
        else:
            return NotImplemented
        return new

    def __neg__(self):
        if self.dtype.kind == 'b':
            raise NotImplementedError('Negating a sparse boolean matrix is not'
                                      ' supported.')
        new = dok_matrix(self.shape, dtype=self.dtype)
        dict.update(new, ((k, -self[k]) for k in iterkeys(self)))
        return new

    def _mul_scalar(self, other):
        res_dtype = upcast_scalar(self.dtype, other)
        # Multiply this scalar by every element.
        new = dok_matrix(self.shape, dtype=res_dtype)
        dict.update(new, ((k, v * other) for k, v in iteritems(self)))
        return new

    def _mul_vector(self, other):
        # matrix * vector
        result = np.zeros(self.shape[0], dtype=upcast(self.dtype, other.dtype))
        for (i, j), v in iteritems(self):
            result[i] += v * other[j]
        return result

    def _mul_multivector(self, other):
        # matrix * multivector
        result_shape = (self.shape[0], other.shape[1])
        result_dtype = upcast(self.dtype, other.dtype)
        result = np.zeros(result_shape, dtype=result_dtype)
        for (i, j), v in iteritems(self):
            result[i,:] += v * other[j,:]
        return result

    def __imul__(self, other):
        if isscalarlike(other):
            dict.update(self, ((k, v * other) for k, v in iteritems(self)))
            return self
        return NotImplemented

    def __truediv__(self, other):
        if isscalarlike(other):
            res_dtype = upcast_scalar(self.dtype, other)
            new = dok_matrix(self.shape, dtype=res_dtype)
            dict.update(new, ((k, v / other) for k, v in iteritems(self)))
            return new
        return self.tocsr() / other

    def __itruediv__(self, other):
        if isscalarlike(other):
            dict.update(self, ((k, v / other) for k, v in iteritems(self)))
            return self
        return NotImplemented

    def __reduce__(self):
        # this approach is necessary because __setstate__ is called after
        # __setitem__ upon unpickling and since __init__ is not called there
        # is no shape attribute hence it is not possible to unpickle it.
        return dict.__reduce__(self)

    # What should len(sparse) return? For consistency with dense matrices,
    # perhaps it should be the number of rows?  For now it returns the number
    # of non-zeros.

    def transpose(self, axes=None, copy=False):
        if axes is not None:
            raise ValueError("Sparse matrices do not support "
                             "an 'axes' parameter because swapping "
                             "dimensions is the only logical permutation.")

        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype, copy=copy)
        dict.update(new, (((right, left), val)
                          for (left, right), val in iteritems(self)))
        return new

    transpose.__doc__ = spmatrix.transpose.__doc__

    def conjtransp(self):
        """Return the conjugate transpose."""
        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype)
        dict.update(new, (((right, left), np.conj(val))
                          for (left, right), val in iteritems(self)))
        return new

    def copy(self):
        new = dok_matrix(self.shape, dtype=self.dtype)
        dict.update(new, self)
        return new

    copy.__doc__ = spmatrix.copy.__doc__

    def getrow(self, i):
        """Returns the i-th row as a (1 x n) DOK matrix."""
        new = dok_matrix((1, self.shape[1]), dtype=self.dtype)
        dict.update(new, (((0, j), self[i, j]) for j in xrange(self.shape[1])))
        return new

    def getcol(self, j):
        """Returns the j-th column as a (m x 1) DOK matrix."""
        new = dok_matrix((self.shape[0], 1), dtype=self.dtype)
        dict.update(new, (((i, 0), self[i, j]) for i in xrange(self.shape[0])))
        return new

    def tocoo(self, copy=False):
        from .coo import coo_matrix
        if self.nnz == 0:
            return coo_matrix(self.shape, dtype=self.dtype)

        idx_dtype = get_index_dtype(maxval=max(self.shape))
        data = np.fromiter(itervalues(self), dtype=self.dtype, count=self.nnz)
        row = np.fromiter((i for i, _ in iterkeys(self)), dtype=idx_dtype, count=self.nnz)
        col = np.fromiter((j for _, j in iterkeys(self)), dtype=idx_dtype, count=self.nnz)
        A = coo_matrix((data, (row, col)), shape=self.shape, dtype=self.dtype)
        A.has_canonical_format = True
        return A

    tocoo.__doc__ = spmatrix.tocoo.__doc__

    def todok(self, copy=False):
        if copy:
            return self.copy()
        return self

    todok.__doc__ = spmatrix.todok.__doc__

    def tocsc(self, copy=False):
        return self.tocoo(copy=False).tocsc(copy=copy)

    tocsc.__doc__ = spmatrix.tocsc.__doc__

    def resize(self, *shape):
        shape = check_shape(shape)
        newM, newN = shape
        M, N = self.shape
        if newM < M or newN < N:
            # Remove all elements outside new dimensions
            for (i, j) in list(iterkeys(self)):
                if i >= newM or j >= newN:
                    del self[i, j]
        self._shape = shape

    resize.__doc__ = spmatrix.resize.__doc__


def isspmatrix_dok(x):
    """Is x of dok_matrix type?

    Parameters
    ----------
    x
        object to check for being a dok matrix

    Returns
    -------
    bool
        True if x is a dok matrix, False otherwise

    Examples
    --------
    >>> from scipy.sparse import dok_matrix, isspmatrix_dok
    >>> isspmatrix_dok(dok_matrix([[5]]))
    True

    >>> from scipy.sparse import dok_matrix, csr_matrix, isspmatrix_dok
    >>> isspmatrix_dok(csr_matrix([[5]]))
    False
    """
    return isinstance(x, dok_matrix)


def _prod(x):
    """Product of a list of numbers; ~40x faster vs np.prod for Python tuples"""
    if len(x) == 0:
        return 1
    return functools.reduce(operator.mul, x)
