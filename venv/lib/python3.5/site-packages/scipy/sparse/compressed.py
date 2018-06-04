"""Base class for sparse matrix formats using compressed storage."""
from __future__ import division, print_function, absolute_import

__all__ = []

from warnings import warn
import operator

import numpy as np
from scipy._lib.six import zip as izip
from scipy._lib._util import _prune_array

from .base import spmatrix, isspmatrix, SparseEfficiencyWarning
from .data import _data_matrix, _minmax_mixin
from .dia import dia_matrix
from . import _sparsetools
from .sputils import (upcast, upcast_char, to_native, isdense, isshape,
                      getdtype, isscalarlike, IndexMixin, get_index_dtype,
                      downcast_intp_index, get_sum_dtype, check_shape)


class _cs_matrix(_data_matrix, _minmax_mixin, IndexMixin):
    """base matrix class for compressed row and column oriented matrices"""

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        _data_matrix.__init__(self)

        if isspmatrix(arg1):
            if arg1.format == self.format and copy:
                arg1 = arg1.copy()
            else:
                arg1 = arg1.asformat(self.format)
            self._set_self(arg1)

        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                # create empty matrix
                self._shape = check_shape(arg1)
                M, N = self.shape
                # Select index dtype large enough to pass array and
                # scalar parameters to sparsetools
                idx_dtype = get_index_dtype(maxval=max(M,N))
                self.data = np.zeros(0, getdtype(dtype, default=float))
                self.indices = np.zeros(0, idx_dtype)
                self.indptr = np.zeros(self._swap((M,N))[0] + 1, dtype=idx_dtype)
            else:
                if len(arg1) == 2:
                    # (data, ij) format
                    from .coo import coo_matrix
                    other = self.__class__(coo_matrix(arg1, shape=shape))
                    self._set_self(other)
                elif len(arg1) == 3:
                    # (data, indices, indptr) format
                    (data, indices, indptr) = arg1

                    # Select index dtype large enough to pass array and
                    # scalar parameters to sparsetools
                    maxval = None
                    if shape is not None:
                        maxval = max(shape)
                    idx_dtype = get_index_dtype((indices, indptr), maxval=maxval, check_contents=True)

                    self.indices = np.array(indices, copy=copy, dtype=idx_dtype)
                    self.indptr = np.array(indptr, copy=copy, dtype=idx_dtype)
                    self.data = np.array(data, copy=copy, dtype=dtype)
                else:
                    raise ValueError("unrecognized %s_matrix constructor usage" %
                            self.format)

        else:
            # must be dense
            try:
                arg1 = np.asarray(arg1)
            except:
                raise ValueError("unrecognized %s_matrix constructor usage" %
                        self.format)
            from .coo import coo_matrix
            self._set_self(self.__class__(coo_matrix(arg1, dtype=dtype)))

        # Read matrix dimensions given, if any
        if shape is not None:
            self._shape = check_shape(shape)
        else:
            if self.shape is None:
                # shape not already set, try to infer dimensions
                try:
                    major_dim = len(self.indptr) - 1
                    minor_dim = self.indices.max() + 1
                except:
                    raise ValueError('unable to infer matrix dimensions')
                else:
                    self._shape = check_shape(self._swap((major_dim,minor_dim)))

        if dtype is not None:
            self.data = np.asarray(self.data, dtype=dtype)

        self.check_format(full_check=False)

    def getnnz(self, axis=None):
        if axis is None:
            return int(self.indptr[-1])
        else:
            if axis < 0:
                axis += 2
            axis, _ = self._swap((axis, 1 - axis))
            _, N = self._swap(self.shape)
            if axis == 0:
                return np.bincount(downcast_intp_index(self.indices),
                                   minlength=N)
            elif axis == 1:
                return np.diff(self.indptr)
            raise ValueError('axis out of bounds')

    getnnz.__doc__ = spmatrix.getnnz.__doc__

    def _set_self(self, other, copy=False):
        """take the member variables of other and assign them to self"""

        if copy:
            other = other.copy()

        self.data = other.data
        self.indices = other.indices
        self.indptr = other.indptr
        self._shape = check_shape(other.shape)

    def check_format(self, full_check=True):
        """check whether the matrix format is valid

        Parameters
        ----------
        full_check : bool, optional
            If `True`, rigorous check, O(N) operations. Otherwise
            basic check, O(1) operations (default True).
        """
        # use _swap to determine proper bounds
        major_name,minor_name = self._swap(('row','column'))
        major_dim,minor_dim = self._swap(self.shape)

        # index arrays should have integer data types
        if self.indptr.dtype.kind != 'i':
            warn("indptr array has non-integer dtype (%s)"
                    % self.indptr.dtype.name)
        if self.indices.dtype.kind != 'i':
            warn("indices array has non-integer dtype (%s)"
                    % self.indices.dtype.name)

        idx_dtype = get_index_dtype((self.indptr, self.indices))
        self.indptr = np.asarray(self.indptr, dtype=idx_dtype)
        self.indices = np.asarray(self.indices, dtype=idx_dtype)
        self.data = to_native(self.data)

        # check array shapes
        if self.data.ndim != 1 or self.indices.ndim != 1 or self.indptr.ndim != 1:
            raise ValueError('data, indices, and indptr should be 1-D')

        # check index pointer
        if (len(self.indptr) != major_dim + 1):
            raise ValueError("index pointer size (%d) should be (%d)" %
                                (len(self.indptr), major_dim + 1))
        if (self.indptr[0] != 0):
            raise ValueError("index pointer should start with 0")

        # check index and data arrays
        if (len(self.indices) != len(self.data)):
            raise ValueError("indices and data should have the same size")
        if (self.indptr[-1] > len(self.indices)):
            raise ValueError("Last value of index pointer should be less than "
                                "the size of index and data arrays")

        self.prune()

        if full_check:
            # check format validity (more expensive)
            if self.nnz > 0:
                if self.indices.max() >= minor_dim:
                    raise ValueError("%s index values must be < %d" %
                                        (minor_name,minor_dim))
                if self.indices.min() < 0:
                    raise ValueError("%s index values must be >= 0" %
                                        minor_name)
                if np.diff(self.indptr).min() < 0:
                    raise ValueError("index pointer values must form a "
                                        "non-decreasing sequence")

        # if not self.has_sorted_indices():
        #    warn('Indices were not in sorted order.  Sorting indices.')
        #    self.sort_indices()
        #    assert(self.has_sorted_indices())
        # TODO check for duplicates?

    #######################
    # Boolean comparisons #
    #######################

    def _scalar_binopt(self, other, op):
        """Scalar version of self._binopt, for cases in which no new nonzeros
        are added. Produces a new spmatrix in canonical form.
        """
        self.sum_duplicates()
        res = self._with_data(op(self.data, other), copy=True)
        res.eliminate_zeros()
        return res

    def __eq__(self, other):
        # Scalar other.
        if isscalarlike(other):
            if np.isnan(other):
                return self.__class__(self.shape, dtype=np.bool_)

            if other == 0:
                warn("Comparing a sparse matrix with 0 using == is inefficient"
                        ", try using != instead.", SparseEfficiencyWarning)
                all_true = self.__class__(np.ones(self.shape, dtype=np.bool_))
                inv = self._scalar_binopt(other, operator.ne)
                return all_true - inv
            else:
                return self._scalar_binopt(other, operator.eq)
        # Dense other.
        elif isdense(other):
            return self.todense() == other
        # Sparse other.
        elif isspmatrix(other):
            warn("Comparing sparse matrices using == is inefficient, try using"
                    " != instead.", SparseEfficiencyWarning)
            #TODO sparse broadcasting
            if self.shape != other.shape:
                return False
            elif self.format != other.format:
                other = other.asformat(self.format)
            res = self._binopt(other,'_ne_')
            all_true = self.__class__(np.ones(self.shape, dtype=np.bool_))
            return all_true - res
        else:
            return False

    def __ne__(self, other):
        # Scalar other.
        if isscalarlike(other):
            if np.isnan(other):
                warn("Comparing a sparse matrix with nan using != is inefficient",
                     SparseEfficiencyWarning)
                all_true = self.__class__(np.ones(self.shape, dtype=np.bool_))
                return all_true
            elif other != 0:
                warn("Comparing a sparse matrix with a nonzero scalar using !="
                     " is inefficient, try using == instead.", SparseEfficiencyWarning)
                all_true = self.__class__(np.ones(self.shape), dtype=np.bool_)
                inv = self._scalar_binopt(other, operator.eq)
                return all_true - inv
            else:
                return self._scalar_binopt(other, operator.ne)
        # Dense other.
        elif isdense(other):
            return self.todense() != other
        # Sparse other.
        elif isspmatrix(other):
            #TODO sparse broadcasting
            if self.shape != other.shape:
                return True
            elif self.format != other.format:
                other = other.asformat(self.format)
            return self._binopt(other,'_ne_')
        else:
            return True

    def _inequality(self, other, op, op_name, bad_scalar_msg):
        # Scalar other.
        if isscalarlike(other):
            if 0 == other and op_name in ('_le_', '_ge_'):
                raise NotImplementedError(" >= and <= don't work with 0.")
            elif op(0, other):
                warn(bad_scalar_msg, SparseEfficiencyWarning)
                other_arr = np.empty(self.shape, dtype=np.result_type(other))
                other_arr.fill(other)
                other_arr = self.__class__(other_arr)
                return self._binopt(other_arr, op_name)
            else:
                return self._scalar_binopt(other, op)
        # Dense other.
        elif isdense(other):
            return op(self.todense(), other)
        # Sparse other.
        elif isspmatrix(other):
            #TODO sparse broadcasting
            if self.shape != other.shape:
                raise ValueError("inconsistent shapes")
            elif self.format != other.format:
                other = other.asformat(self.format)
            if op_name not in ('_ge_', '_le_'):
                return self._binopt(other, op_name)

            warn("Comparing sparse matrices using >= and <= is inefficient, "
                 "using <, >, or !=, instead.", SparseEfficiencyWarning)
            all_true = self.__class__(np.ones(self.shape))
            res = self._binopt(other, '_gt_' if op_name == '_le_' else '_lt_')
            return all_true - res
        else:
            raise ValueError("Operands could not be compared.")

    def __lt__(self, other):
        return self._inequality(other, operator.lt, '_lt_',
                                "Comparing a sparse matrix with a scalar "
                                "greater than zero using < is inefficient, "
                                "try using >= instead.")

    def __gt__(self, other):
        return self._inequality(other, operator.gt, '_gt_',
                                "Comparing a sparse matrix with a scalar "
                                "less than zero using > is inefficient, "
                                "try using <= instead.")

    def __le__(self, other):
        return self._inequality(other, operator.le, '_le_',
                                "Comparing a sparse matrix with a scalar "
                                "greater than zero using <= is inefficient, "
                                "try using > instead.")

    def __ge__(self,other):
        return self._inequality(other, operator.ge, '_ge_',
                                "Comparing a sparse matrix with a scalar "
                                "less than zero using >= is inefficient, "
                                "try using < instead.")

    #################################
    # Arithmetic operator overrides #
    #################################

    def _add_dense(self, other):
        if other.shape != self.shape:
            raise ValueError('Incompatible shapes.')
        dtype = upcast_char(self.dtype.char, other.dtype.char)
        order = self._swap('CF')[0]
        result = np.array(other, dtype=dtype, order=order, copy=True)
        M, N = self._swap(self.shape)
        y = result if result.flags.c_contiguous else result.T
        _sparsetools.csr_todense(M, N, self.indptr, self.indices, self.data, y)
        return np.matrix(result, copy=False)

    def _add_sparse(self, other):
        return self._binopt(other, '_plus_')

    def _sub_sparse(self, other):
        return self._binopt(other, '_minus_')

    def multiply(self, other):
        """Point-wise multiplication by another matrix, vector, or
        scalar.
        """
        # Scalar multiplication.
        if isscalarlike(other):
            return self._mul_scalar(other)
        # Sparse matrix or vector.
        if isspmatrix(other):
            if self.shape == other.shape:
                other = self.__class__(other)
                return self._binopt(other, '_elmul_')
            # Single element.
            elif other.shape == (1,1):
                return self._mul_scalar(other.toarray()[0, 0])
            elif self.shape == (1,1):
                return other._mul_scalar(self.toarray()[0, 0])
            # A row times a column.
            elif self.shape[1] == 1 and other.shape[0] == 1:
                return self._mul_sparse_matrix(other.tocsc())
            elif self.shape[0] == 1 and other.shape[1] == 1:
                return other._mul_sparse_matrix(self.tocsc())
            # Row vector times matrix. other is a row.
            elif other.shape[0] == 1 and self.shape[1] == other.shape[1]:
                other = dia_matrix((other.toarray().ravel(), [0]),
                                   shape=(other.shape[1], other.shape[1]))
                return self._mul_sparse_matrix(other)
            # self is a row.
            elif self.shape[0] == 1 and self.shape[1] == other.shape[1]:
                copy = dia_matrix((self.toarray().ravel(), [0]),
                                  shape=(self.shape[1], self.shape[1]))
                return other._mul_sparse_matrix(copy)
            # Column vector times matrix. other is a column.
            elif other.shape[1] == 1 and self.shape[0] == other.shape[0]:
                other = dia_matrix((other.toarray().ravel(), [0]),
                                   shape=(other.shape[0], other.shape[0]))
                return other._mul_sparse_matrix(self)
            # self is a column.
            elif self.shape[1] == 1 and self.shape[0] == other.shape[0]:
                copy = dia_matrix((self.toarray().ravel(), [0]),
                                  shape=(self.shape[0], self.shape[0]))
                return copy._mul_sparse_matrix(other)
            else:
                raise ValueError("inconsistent shapes")

        # Assume other is a dense matrix/array, which produces a single-item
        # object array if other isn't convertible to ndarray.
        other = np.atleast_2d(other)

        if other.ndim != 2:
            return np.multiply(self.toarray(), other)
        # Single element / wrapped object.
        if other.size == 1:
            return self._mul_scalar(other.flat[0])
        # Fast case for trivial sparse matrix.
        elif self.shape == (1, 1):
            return np.multiply(self.toarray()[0,0], other)

        from .coo import coo_matrix
        ret = self.tocoo()
        # Matching shapes.
        if self.shape == other.shape:
            data = np.multiply(ret.data, other[ret.row, ret.col])
        # Sparse row vector times...
        elif self.shape[0] == 1:
            if other.shape[1] == 1:  # Dense column vector.
                data = np.multiply(ret.data, other)
            elif other.shape[1] == self.shape[1]:  # Dense matrix.
                data = np.multiply(ret.data, other[:, ret.col])
            else:
                raise ValueError("inconsistent shapes")
            row = np.repeat(np.arange(other.shape[0]), len(ret.row))
            col = np.tile(ret.col, other.shape[0])
            return coo_matrix((data.view(np.ndarray).ravel(), (row, col)),
                              shape=(other.shape[0], self.shape[1]),
                              copy=False)
        # Sparse column vector times...
        elif self.shape[1] == 1:
            if other.shape[0] == 1:  # Dense row vector.
                data = np.multiply(ret.data[:, None], other)
            elif other.shape[0] == self.shape[0]:  # Dense matrix.
                data = np.multiply(ret.data[:, None], other[ret.row])
            else:
                raise ValueError("inconsistent shapes")
            row = np.repeat(ret.row, other.shape[1])
            col = np.tile(np.arange(other.shape[1]), len(ret.col))
            return coo_matrix((data.view(np.ndarray).ravel(), (row, col)),
                              shape=(self.shape[0], other.shape[1]),
                              copy=False)
        # Sparse matrix times dense row vector.
        elif other.shape[0] == 1 and self.shape[1] == other.shape[1]:
            data = np.multiply(ret.data, other[:, ret.col].ravel())
        # Sparse matrix times dense column vector.
        elif other.shape[1] == 1 and self.shape[0] == other.shape[0]:
            data = np.multiply(ret.data, other[ret.row].ravel())
        else:
            raise ValueError("inconsistent shapes")
        ret.data = data.view(np.ndarray).ravel()
        return ret

    ###########################
    # Multiplication handlers #
    ###########################

    def _mul_vector(self, other):
        M,N = self.shape

        # output array
        result = np.zeros(M, dtype=upcast_char(self.dtype.char,
                                               other.dtype.char))

        # csr_matvec or csc_matvec
        fn = getattr(_sparsetools,self.format + '_matvec')
        fn(M, N, self.indptr, self.indices, self.data, other, result)

        return result

    def _mul_multivector(self, other):
        M,N = self.shape
        n_vecs = other.shape[1]  # number of column vectors

        result = np.zeros((M,n_vecs), dtype=upcast_char(self.dtype.char,
                                                        other.dtype.char))

        # csr_matvecs or csc_matvecs
        fn = getattr(_sparsetools,self.format + '_matvecs')
        fn(M, N, n_vecs, self.indptr, self.indices, self.data, other.ravel(), result.ravel())

        return result

    def _mul_sparse_matrix(self, other):
        M, K1 = self.shape
        K2, N = other.shape

        major_axis = self._swap((M,N))[0]
        other = self.__class__(other)  # convert to this format

        idx_dtype = get_index_dtype((self.indptr, self.indices,
                                     other.indptr, other.indices),
                                    maxval=M*N)
        indptr = np.empty(major_axis + 1, dtype=idx_dtype)

        fn = getattr(_sparsetools, self.format + '_matmat_pass1')
        fn(M, N,
           np.asarray(self.indptr, dtype=idx_dtype),
           np.asarray(self.indices, dtype=idx_dtype),
           np.asarray(other.indptr, dtype=idx_dtype),
           np.asarray(other.indices, dtype=idx_dtype),
           indptr)

        nnz = indptr[-1]
        idx_dtype = get_index_dtype((self.indptr, self.indices,
                                     other.indptr, other.indices),
                                    maxval=nnz)
        indptr = np.asarray(indptr, dtype=idx_dtype)
        indices = np.empty(nnz, dtype=idx_dtype)
        data = np.empty(nnz, dtype=upcast(self.dtype, other.dtype))

        fn = getattr(_sparsetools, self.format + '_matmat_pass2')
        fn(M, N, np.asarray(self.indptr, dtype=idx_dtype),
           np.asarray(self.indices, dtype=idx_dtype),
           self.data,
           np.asarray(other.indptr, dtype=idx_dtype),
           np.asarray(other.indices, dtype=idx_dtype),
           other.data,
           indptr, indices, data)

        return self.__class__((data,indices,indptr),shape=(M,N))

    def diagonal(self, k=0):
        rows, cols = self.shape
        if k <= -rows or k >= cols:
            raise ValueError("k exceeds matrix dimensions")
        fn = getattr(_sparsetools, self.format + "_diagonal")
        y = np.empty(min(rows + min(k, 0), cols - max(k, 0)),
                     dtype=upcast(self.dtype))
        fn(k, self.shape[0], self.shape[1], self.indptr, self.indices,
           self.data, y)
        return y

    diagonal.__doc__ = spmatrix.diagonal.__doc__

    #####################
    # Other binary ops  #
    #####################

    def _maximum_minimum(self, other, npop, op_name, dense_check):
        if isscalarlike(other):
            if dense_check(other):
                warn("Taking maximum (minimum) with > 0 (< 0) number results to "
                     "a dense matrix.",
                     SparseEfficiencyWarning)
                other_arr = np.empty(self.shape, dtype=np.asarray(other).dtype)
                other_arr.fill(other)
                other_arr = self.__class__(other_arr)
                return self._binopt(other_arr, op_name)
            else:
                self.sum_duplicates()
                new_data = npop(self.data, np.asarray(other))
                mat = self.__class__((new_data, self.indices, self.indptr),
                                     dtype=new_data.dtype, shape=self.shape)
                return mat
        elif isdense(other):
            return npop(self.todense(), other)
        elif isspmatrix(other):
            return self._binopt(other, op_name)
        else:
            raise ValueError("Operands not compatible.")

    def maximum(self, other):
        return self._maximum_minimum(other, np.maximum, '_maximum_', lambda x: np.asarray(x) > 0)

    maximum.__doc__ = spmatrix.maximum.__doc__

    def minimum(self, other):
        return self._maximum_minimum(other, np.minimum, '_minimum_', lambda x: np.asarray(x) < 0)

    minimum.__doc__ = spmatrix.minimum.__doc__

    #####################
    # Reduce operations #
    #####################

    def sum(self, axis=None, dtype=None, out=None):
        """Sum the matrix over the given axis.  If the axis is None, sum
        over both rows and columns, returning a scalar.
        """
        # The spmatrix base class already does axis=0 and axis=1 efficiently
        # so we only do the case axis=None here
        if (not hasattr(self, 'blocksize') and
                    axis in self._swap(((1, -1), (0, 2)))[0]):
            # faster than multiplication for large minor axis in CSC/CSR
            res_dtype = get_sum_dtype(self.dtype)
            ret = np.zeros(len(self.indptr) - 1, dtype=res_dtype)

            major_index, value = self._minor_reduce(np.add)
            ret[major_index] = value
            ret = np.asmatrix(ret)
            if axis % 2 == 1:
                ret = ret.T

            if out is not None and out.shape != ret.shape:
                raise ValueError('dimensions do not match')

            return ret.sum(axis=(), dtype=dtype, out=out)
        # spmatrix will handle the remaining situations when axis
        # is in {None, -1, 0, 1}
        else:
            return spmatrix.sum(self, axis=axis, dtype=dtype, out=out)

    sum.__doc__ = spmatrix.sum.__doc__

    def _minor_reduce(self, ufunc, data=None):
        """Reduce nonzeros with a ufunc over the minor axis when non-empty

        Can be applied to a function of self.data by supplying data parameter.

        Warning: this does not call sum_duplicates()

        Returns
        -------
        major_index : array of ints
            Major indices where nonzero

        value : array of self.dtype
            Reduce result for nonzeros in each major_index
        """
        if data is None:
            data = self.data
        major_index = np.flatnonzero(np.diff(self.indptr))
        value = ufunc.reduceat(data,
                               downcast_intp_index(self.indptr[major_index]))
        return major_index, value

    #######################
    # Getting and Setting #
    #######################

    def __setitem__(self, index, x):
        # Process arrays from IndexMixin
        i, j = self._unpack_index(index)
        i, j = self._index_to_arrays(i, j)

        if isspmatrix(x):
            broadcast_row = x.shape[0] == 1 and i.shape[0] != 1
            broadcast_col = x.shape[1] == 1 and i.shape[1] != 1
            if not ((broadcast_row or x.shape[0] == i.shape[0]) and
                    (broadcast_col or x.shape[1] == i.shape[1])):
                raise ValueError("shape mismatch in assignment")

            # clear entries that will be overwritten
            ci, cj = self._swap((i.ravel(), j.ravel()))
            self._zero_many(ci, cj)

            x = x.tocoo(copy=True)
            x.sum_duplicates()
            r, c = x.row, x.col
            x = np.asarray(x.data, dtype=self.dtype)
            if broadcast_row:
                r = np.repeat(np.arange(i.shape[0]), len(r))
                c = np.tile(c, i.shape[0])
                x = np.tile(x, i.shape[0])
            if broadcast_col:
                r = np.repeat(r, i.shape[1])
                c = np.tile(np.arange(i.shape[1]), len(c))
                x = np.repeat(x, i.shape[1])
            # only assign entries in the new sparsity structure
            i = i[r, c]
            j = j[r, c]
        else:
            # Make x and i into the same shape
            x = np.asarray(x, dtype=self.dtype)
            x, _ = np.broadcast_arrays(x, i)

            if x.shape != i.shape:
                raise ValueError("shape mismatch in assignment")

        if np.size(x) == 0:
            return
        i, j = self._swap((i.ravel(), j.ravel()))
        self._set_many(i, j, x.ravel())

    def _setdiag(self, values, k):
        if 0 in self.shape:
            return

        M, N = self.shape
        broadcast = (values.ndim == 0)

        if k < 0:
            if broadcast:
                max_index = min(M + k, N)
            else:
                max_index = min(M + k, N, len(values))
            i = np.arange(max_index, dtype=self.indices.dtype)
            j = np.arange(max_index, dtype=self.indices.dtype)
            i -= k

        else:
            if broadcast:
                max_index = min(M, N - k)
            else:
                max_index = min(M, N - k, len(values))
            i = np.arange(max_index, dtype=self.indices.dtype)
            j = np.arange(max_index, dtype=self.indices.dtype)
            j += k

        if not broadcast:
            values = values[:len(i)]

        self[i, j] = values

    def _prepare_indices(self, i, j):
        M, N = self._swap(self.shape)

        def check_bounds(indices, bound):
            idx = indices.max()
            if idx >= bound:
                raise IndexError('index (%d) out of range (>= %d)' %
                                 (idx, bound))
            idx = indices.min()
            if idx < -bound:
                raise IndexError('index (%d) out of range (< -%d)' %
                                 (idx, bound))

        check_bounds(i, M)
        check_bounds(j, N)

        i = np.asarray(i, dtype=self.indices.dtype)
        j = np.asarray(j, dtype=self.indices.dtype)
        return i, j, M, N

    def _set_many(self, i, j, x):
        """Sets value at each (i, j) to x

        Here (i,j) index major and minor respectively, and must not contain
        duplicate entries.
        """
        i, j, M, N = self._prepare_indices(i, j)

        n_samples = len(x)
        offsets = np.empty(n_samples, dtype=self.indices.dtype)
        ret = _sparsetools.csr_sample_offsets(M, N, self.indptr, self.indices,
                                              n_samples, i, j, offsets)
        if ret == 1:
            # rinse and repeat
            self.sum_duplicates()
            _sparsetools.csr_sample_offsets(M, N, self.indptr,
                                            self.indices, n_samples, i, j,
                                            offsets)

        if -1 not in offsets:
            # only affects existing non-zero cells
            self.data[offsets] = x
            return

        else:
            warn("Changing the sparsity structure of a %s_matrix is expensive. "
                 "lil_matrix is more efficient." % self.format,
                 SparseEfficiencyWarning)
            # replace where possible
            mask = offsets > -1
            self.data[offsets[mask]] = x[mask]
            # only insertions remain
            mask = ~mask
            i = i[mask]
            i[i < 0] += M
            j = j[mask]
            j[j < 0] += N
            self._insert_many(i, j, x[mask])

    def _zero_many(self, i, j):
        """Sets value at each (i, j) to zero, preserving sparsity structure.

        Here (i,j) index major and minor respectively.
        """
        i, j, M, N = self._prepare_indices(i, j)

        n_samples = len(i)
        offsets = np.empty(n_samples, dtype=self.indices.dtype)
        ret = _sparsetools.csr_sample_offsets(M, N, self.indptr, self.indices,
                                              n_samples, i, j, offsets)
        if ret == 1:
            # rinse and repeat
            self.sum_duplicates()
            _sparsetools.csr_sample_offsets(M, N, self.indptr,
                                            self.indices, n_samples, i, j,
                                            offsets)

        # only assign zeros to the existing sparsity structure
        self.data[offsets[offsets > -1]] = 0

    def _insert_many(self, i, j, x):
        """Inserts new nonzero at each (i, j) with value x

        Here (i,j) index major and minor respectively.
        i, j and x must be non-empty, 1d arrays.
        Inserts each major group (e.g. all entries per row) at a time.
        Maintains has_sorted_indices property.
        Modifies i, j, x in place.
        """
        order = np.argsort(i, kind='mergesort')  # stable for duplicates
        i = i.take(order, mode='clip')
        j = j.take(order, mode='clip')
        x = x.take(order, mode='clip')

        do_sort = self.has_sorted_indices

        # Update index data type
        idx_dtype = get_index_dtype((self.indices, self.indptr),
                                    maxval=(self.indptr[-1] + x.size))
        self.indptr = np.asarray(self.indptr, dtype=idx_dtype)
        self.indices = np.asarray(self.indices, dtype=idx_dtype)
        i = np.asarray(i, dtype=idx_dtype)
        j = np.asarray(j, dtype=idx_dtype)

        # Collate old and new in chunks by major index
        indices_parts = []
        data_parts = []
        ui, ui_indptr = np.unique(i, return_index=True)
        ui_indptr = np.append(ui_indptr, len(j))
        new_nnzs = np.diff(ui_indptr)
        prev = 0
        for c, (ii, js, je) in enumerate(izip(ui, ui_indptr, ui_indptr[1:])):
            # old entries
            start = self.indptr[prev]
            stop = self.indptr[ii]
            indices_parts.append(self.indices[start:stop])
            data_parts.append(self.data[start:stop])

            # handle duplicate j: keep last setting
            uj, uj_indptr = np.unique(j[js:je][::-1], return_index=True)
            if len(uj) == je - js:
                indices_parts.append(j[js:je])
                data_parts.append(x[js:je])
            else:
                indices_parts.append(j[js:je][::-1][uj_indptr])
                data_parts.append(x[js:je][::-1][uj_indptr])
                new_nnzs[c] = len(uj)

            prev = ii

        # remaining old entries
        start = self.indptr[ii]
        indices_parts.append(self.indices[start:])
        data_parts.append(self.data[start:])

        # update attributes
        self.indices = np.concatenate(indices_parts)
        self.data = np.concatenate(data_parts)
        nnzs = np.asarray(np.ediff1d(self.indptr, to_begin=0), dtype=idx_dtype)
        nnzs[1:][ui] += new_nnzs
        self.indptr = np.cumsum(nnzs, out=nnzs)

        if do_sort:
            # TODO: only sort where necessary
            self.has_sorted_indices = False
            self.sort_indices()

        self.check_format(full_check=False)

    def _get_single_element(self, row, col):
        M, N = self.shape
        if (row < 0):
            row += M
        if (col < 0):
            col += N
        if not (0 <= row < M) or not (0 <= col < N):
            raise IndexError("index out of bounds: 0<=%d<%d, 0<=%d<%d" %
                             (row, M, col, N))

        major_index, minor_index = self._swap((row, col))

        start = self.indptr[major_index]
        end = self.indptr[major_index + 1]

        if self.has_sorted_indices:
            # Copies may be made, if dtypes of indices are not identical
            minor_index = self.indices.dtype.type(minor_index)
            minor_indices = self.indices[start:end]
            insert_pos_left = np.searchsorted(
                minor_indices, minor_index, side='left')
            insert_pos_right = insert_pos_left + np.searchsorted(
                minor_indices[insert_pos_left:], minor_index, side='right')
            return self.data[start + insert_pos_left:
                             start + insert_pos_right].sum(dtype=self.dtype)
        else:
            return np.compress(minor_index == self.indices[start:end],
                               self.data[start:end]).sum(dtype=self.dtype)

    def _get_submatrix(self, slice0, slice1):
        """Return a submatrix of this matrix (new matrix is created)."""

        slice0, slice1 = self._swap((slice0,slice1))
        shape0, shape1 = self._swap(self.shape)

        def _process_slice(sl, num):
            if isinstance(sl, slice):
                i0, i1 = sl.start, sl.stop
                if i0 is None:
                    i0 = 0
                elif i0 < 0:
                    i0 = num + i0

                if i1 is None:
                    i1 = num
                elif i1 < 0:
                    i1 = num + i1

                return i0, i1

            elif np.isscalar(sl):
                if sl < 0:
                    sl += num

                return sl, sl + 1

            else:
                return sl[0], sl[1]

        def _in_bounds(i0, i1, num):
            if not (0 <= i0 < num) or not (0 < i1 <= num) or not (i0 < i1):
                raise IndexError("index out of bounds: 0<=%d<%d, 0<=%d<%d, %d<%d" %
                                    (i0, num, i1, num, i0, i1))

        i0, i1 = _process_slice(slice0, shape0)
        j0, j1 = _process_slice(slice1, shape1)
        _in_bounds(i0, i1, shape0)
        _in_bounds(j0, j1, shape1)

        aux = _sparsetools.get_csr_submatrix(shape0, shape1,
                                             self.indptr, self.indices,
                                             self.data,
                                             i0, i1, j0, j1)

        data, indices, indptr = aux[2], aux[1], aux[0]
        shape = self._swap((i1 - i0, j1 - j0))

        return self.__class__((data, indices, indptr), shape=shape)

    ######################
    # Conversion methods #
    ######################

    def tocoo(self, copy=True):
        major_dim, minor_dim = self._swap(self.shape)
        minor_indices = self.indices
        major_indices = np.empty(len(minor_indices), dtype=self.indices.dtype)
        _sparsetools.expandptr(major_dim, self.indptr, major_indices)
        row, col = self._swap((major_indices, minor_indices))

        from .coo import coo_matrix
        return coo_matrix((self.data, (row, col)), self.shape, copy=copy,
                          dtype=self.dtype)

    tocoo.__doc__ = spmatrix.tocoo.__doc__

    def toarray(self, order=None, out=None):
        if out is None and order is None:
            order = self._swap('cf')[0]
        out = self._process_toarray_args(order, out)
        if not (out.flags.c_contiguous or out.flags.f_contiguous):
            raise ValueError('Output array must be C or F contiguous')
        # align ideal order with output array order
        if out.flags.c_contiguous:
            x = self.tocsr()
            y = out
        else:
            x = self.tocsc()
            y = out.T
        M, N = x._swap(x.shape)
        _sparsetools.csr_todense(M, N, x.indptr, x.indices, x.data, y)
        return out

    toarray.__doc__ = spmatrix.toarray.__doc__

    ##############################################################
    # methods that examine or modify the internal data structure #
    ##############################################################

    def eliminate_zeros(self):
        """Remove zero entries from the matrix

        This is an *in place* operation
        """
        M, N = self._swap(self.shape)
        _sparsetools.csr_eliminate_zeros(M, N, self.indptr, self.indices,
                                         self.data)
        self.prune()  # nnz may have changed

    def __get_has_canonical_format(self):
        """Determine whether the matrix has sorted indices and no duplicates

        Returns
            - True: if the above applies
            - False: otherwise

        has_canonical_format implies has_sorted_indices, so if the latter flag
        is False, so will the former be; if the former is found True, the
        latter flag is also set.
        """

        # first check to see if result was cached
        if not getattr(self, '_has_sorted_indices', True):
            # not sorted => not canonical
            self._has_canonical_format = False
        elif not hasattr(self, '_has_canonical_format'):
            self.has_canonical_format = _sparsetools.csr_has_canonical_format(
                len(self.indptr) - 1, self.indptr, self.indices)
        return self._has_canonical_format

    def __set_has_canonical_format(self, val):
        self._has_canonical_format = bool(val)
        if val:
            self.has_sorted_indices = True

    has_canonical_format = property(fget=__get_has_canonical_format,
                                    fset=__set_has_canonical_format)

    def sum_duplicates(self):
        """Eliminate duplicate matrix entries by adding them together

        The is an *in place* operation
        """
        if self.has_canonical_format:
            return
        self.sort_indices()

        M, N = self._swap(self.shape)
        _sparsetools.csr_sum_duplicates(M, N, self.indptr, self.indices,
                                        self.data)

        self.prune()  # nnz may have changed
        self.has_canonical_format = True

    def __get_sorted(self):
        """Determine whether the matrix has sorted indices

        Returns
            - True: if the indices of the matrix are in sorted order
            - False: otherwise

        """

        # first check to see if result was cached
        if not hasattr(self,'_has_sorted_indices'):
            self._has_sorted_indices = _sparsetools.csr_has_sorted_indices(
                len(self.indptr) - 1, self.indptr, self.indices)
        return self._has_sorted_indices

    def __set_sorted(self, val):
        self._has_sorted_indices = bool(val)

    has_sorted_indices = property(fget=__get_sorted, fset=__set_sorted)

    def sorted_indices(self):
        """Return a copy of this matrix with sorted indices
        """
        A = self.copy()
        A.sort_indices()
        return A

        # an alternative that has linear complexity is the following
        # although the previous option is typically faster
        # return self.toother().toother()

    def sort_indices(self):
        """Sort the indices of this matrix *in place*
        """

        if not self.has_sorted_indices:
            _sparsetools.csr_sort_indices(len(self.indptr) - 1, self.indptr,
                                          self.indices, self.data)
            self.has_sorted_indices = True

    def prune(self):
        """Remove empty space after all non-zero elements.
        """
        major_dim = self._swap(self.shape)[0]

        if len(self.indptr) != major_dim + 1:
            raise ValueError('index pointer has invalid length')
        if len(self.indices) < self.nnz:
            raise ValueError('indices array has fewer than nnz elements')
        if len(self.data) < self.nnz:
            raise ValueError('data array has fewer than nnz elements')

        self.indices = _prune_array(self.indices[:self.nnz])
        self.data = _prune_array(self.data[:self.nnz])

    def resize(self, *shape):
        shape = check_shape(shape)
        if hasattr(self, 'blocksize'):
            bm, bn = self.blocksize
            new_M, rm = divmod(shape[0], bm)
            new_N, rn = divmod(shape[1], bn)
            if rm or rn:
                raise ValueError("shape must be divisible into %s blocks. "
                                 "Got %s" % (self.blocksize, shape))
            M, N = self.shape[0] // bm, self.shape[1] // bn
        else:
            new_M, new_N = self._swap(shape)
            M, N = self._swap(self.shape)

        if new_M < M:
            self.indices = self.indices[:self.indptr[new_M]]
            self.data = self.data[:self.indptr[new_M]]
            self.indptr = self.indptr[:new_M + 1]
        elif new_M > M:
            self.indptr = np.resize(self.indptr, new_M + 1)
            self.indptr[M + 1:].fill(self.indptr[M])

        if new_N < N:
            mask = self.indices < new_N
            if not np.all(mask):
                self.indices = self.indices[mask]
                self.data = self.data[mask]
                major_index, val = self._minor_reduce(np.add, mask)
                self.indptr.fill(0)
                self.indptr[1:][major_index] = val
                np.cumsum(self.indptr, out=self.indptr)

        self._shape = shape

    resize.__doc__ = spmatrix.resize.__doc__

    ###################
    # utility methods #
    ###################

    # needed by _data_matrix
    def _with_data(self,data,copy=True):
        """Returns a matrix with the same sparsity structure as self,
        but with different data.  By default the structure arrays
        (i.e. .indptr and .indices) are copied.
        """
        if copy:
            return self.__class__((data,self.indices.copy(),self.indptr.copy()),
                                   shape=self.shape,dtype=data.dtype)
        else:
            return self.__class__((data,self.indices,self.indptr),
                                   shape=self.shape,dtype=data.dtype)

    def _binopt(self, other, op):
        """apply the binary operation fn to two sparse matrices."""
        other = self.__class__(other)

        # e.g. csr_plus_csr, csr_minus_csr, etc.
        fn = getattr(_sparsetools, self.format + op + self.format)

        maxnnz = self.nnz + other.nnz
        idx_dtype = get_index_dtype((self.indptr, self.indices,
                                     other.indptr, other.indices),
                                    maxval=maxnnz)
        indptr = np.empty(self.indptr.shape, dtype=idx_dtype)
        indices = np.empty(maxnnz, dtype=idx_dtype)

        bool_ops = ['_ne_', '_lt_', '_gt_', '_le_', '_ge_']
        if op in bool_ops:
            data = np.empty(maxnnz, dtype=np.bool_)
        else:
            data = np.empty(maxnnz, dtype=upcast(self.dtype, other.dtype))

        fn(self.shape[0], self.shape[1],
           np.asarray(self.indptr, dtype=idx_dtype),
           np.asarray(self.indices, dtype=idx_dtype),
           self.data,
           np.asarray(other.indptr, dtype=idx_dtype),
           np.asarray(other.indices, dtype=idx_dtype),
           other.data,
           indptr, indices, data)

        A = self.__class__((data, indices, indptr), shape=self.shape)
        A.prune()

        return A

    def _divide_sparse(self, other):
        """
        Divide this matrix by a second sparse matrix.
        """
        if other.shape != self.shape:
            raise ValueError('inconsistent shapes')

        r = self._binopt(other, '_eldiv_')

        if np.issubdtype(r.dtype, np.inexact):
            # Eldiv leaves entries outside the combined sparsity
            # pattern empty, so they must be filled manually.
            # Everything outside of other's sparsity is NaN, and everything
            # inside it is either zero or defined by eldiv.
            out = np.empty(self.shape, dtype=self.dtype)
            out.fill(np.nan)
            row, col = other.nonzero()
            out[row, col] = 0
            r = r.tocoo()
            out[r.row, r.col] = r.data
            out = np.matrix(out)
        else:
            # integers types go with nan <-> 0
            out = r

        return out
