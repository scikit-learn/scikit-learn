# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Lars Buitinck
#         Danny Sullivan <dsullivan7@hotmail.com>
#
# License: BSD 3 clause

cimport cython
from libc.limits cimport INT_MAX
from libc.math cimport sqrt
import numpy as np
cimport numpy as np

from ._cython_blas cimport _dot, _scal, _axpy


np.import_array()


cdef class WeightVector(object):
    """Dense vector represented by a scalar and a numpy array.

    The class provides methods to ``add`` a sparse vector
    and scale the vector.
    Representing a vector explicitly as a scalar times a
    vector allows for efficient scaling operations.

    Attributes
    ----------
    w : ndarray, dtype=double, order='C'
        The numpy array which backs the weight vector.
    aw : ndarray, dtype=double, order='C'
        The numpy array which backs the average_weight vector.
    w_data_ptr : double*
        A pointer to the data of the numpy array.
    wscale : double
        The scale of the vector.
    n_features : int
        The number of features (= dimensionality of ``w``).
    sq_norm : double
        The squared norm of ``w``.
    """

    def __cinit__(self,
                  np.ndarray[double, ndim=1, mode='c'] w,
                  np.ndarray[double, ndim=1, mode='c'] aw):
        cdef double *wdata = <double *>w.data

        if w.shape[0] > INT_MAX:
            raise ValueError("More than %d features not supported; got %d."
                             % (INT_MAX, w.shape[0]))
        self.w = w
        self.w_data_ptr = wdata
        self.wscale = 1.0
        self.n_features = w.shape[0]
        self.sq_norm = _dot(<int>w.shape[0], wdata, 1, wdata, 1)

        self.aw = aw
        if self.aw is not None:
            self.aw_data_ptr = <double *>aw.data
            self.average_a = 0.0
            self.average_b = 1.0

    cdef void add(self, double *x_data_ptr, int *x_ind_ptr, int xnnz,
                  double c) nogil:
        """Scales sample x by constant c and adds it to the weight vector.

        This operation updates ``sq_norm``.

        Parameters
        ----------
        x_data_ptr : double*
            The array which holds the feature values of ``x``.
        x_ind_ptr : np.intc*
            The array which holds the feature indices of ``x``.
        xnnz : int
            The number of non-zero features of ``x``.
        c : double
            The scaling constant for the example.
        """
        cdef int j
        cdef int idx
        cdef double val
        cdef double innerprod = 0.0
        cdef double xsqnorm = 0.0

        # the next two lines save a factor of 2!
        cdef double wscale = self.wscale
        cdef double* w_data_ptr = self.w_data_ptr

        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            innerprod += (w_data_ptr[idx] * val)
            xsqnorm += (val * val)
            w_data_ptr[idx] += val * (c / wscale)

        self.sq_norm += (xsqnorm * c * c) + (2.0 * innerprod * wscale * c)

    # Update the average weights according to the sparse trick defined
    # here: https://research.microsoft.com/pubs/192769/tricks-2012.pdf
    # by Leon Bottou
    cdef void add_average(self, double *x_data_ptr, int *x_ind_ptr, int xnnz,
                          double c, double num_iter) nogil:
        """Updates the average weight vector.

        Parameters
        ----------
        x_data_ptr : double*
            The array which holds the feature values of ``x``.
        x_ind_ptr : np.intc*
            The array which holds the feature indices of ``x``.
        xnnz : int
            The number of non-zero features of ``x``.
        c : double
            The scaling constant for the example.
        num_iter : double
            The total number of iterations.
        """
        cdef int j
        cdef int idx
        cdef double val
        cdef double mu = 1.0 / num_iter
        cdef double average_a = self.average_a
        cdef double wscale = self.wscale
        cdef double* aw_data_ptr = self.aw_data_ptr

        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            aw_data_ptr[idx] += (self.average_a * val * (-c / wscale))

        # Once the sample has been processed
        # update the average_a and average_b
        if num_iter > 1:
            self.average_b /= (1.0 - mu)
        self.average_a += mu * self.average_b * wscale

    cdef double dot(self, double *x_data_ptr, int *x_ind_ptr,
                    int xnnz) nogil:
        """Computes the dot product of a sample x and the weight vector.

        Parameters
        ----------
        x_data_ptr : double*
            The array which holds the feature values of ``x``.
        x_ind_ptr : np.intc*
            The array which holds the feature indices of ``x``.
        xnnz : int
            The number of non-zero features of ``x`` (length of x_ind_ptr).

        Returns
        -------
        innerprod : double
            The inner product of ``x`` and ``w``.
        """
        cdef int j
        cdef int idx
        cdef double innerprod = 0.0
        cdef double* w_data_ptr = self.w_data_ptr
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            innerprod += w_data_ptr[idx] * x_data_ptr[j]
        innerprod *= self.wscale
        return innerprod

    cdef void scale(self, double c) nogil:
        """Scales the weight vector by a constant ``c``.

        It updates ``wscale`` and ``sq_norm``. If ``wscale`` gets too
        small we call ``reset_swcale``."""
        self.wscale *= c
        self.sq_norm *= (c * c)
        if self.wscale < 1e-9:
            self.reset_wscale()

    cdef void reset_wscale(self) nogil:
        """Scales each coef of ``w`` by ``wscale`` and resets it to 1. """
        if self.aw is not None:
            _axpy(<int>self.aw.shape[0], self.average_a,
                  <double *>self.w.data, 1, <double *>self.aw.data, 1)
            _scal(<int>self.aw.shape[0], 1.0 / self.average_b,
                  <double *>self.aw.data, 1)
            self.average_a = 0.0
            self.average_b = 1.0

        _scal(<int>self.w.shape[0], self.wscale, <double *>self.w.data, 1)
        self.wscale = 1.0

    cdef double norm(self) nogil:
        """The L2 norm of the weight vector. """
        return sqrt(self.sq_norm)
