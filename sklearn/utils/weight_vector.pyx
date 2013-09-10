# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# Licence: BSD 3 clause

cimport cython
from libc.limits cimport INT_MAX
from libc.math cimport sqrt
import numpy as np
cimport numpy as np

np.import_array()

cdef extern from "cblas.h":
    void dscal "cblas_dscal"(int, double, double *, int) nogil


cdef class WeightVector(object):
    """Dense vector represented by a scalar and a numpy array.

    The class provides methods to ``add`` a sparse vector
    and scale the vector.
    Representing a vector explicitly as a scalar times a
    vector allows for efficient scaling operations.

    Attributes
    ----------
    w : ndarray, dtype=np.float64, order='C'
        The numpy array which backs the weight vector.
    w_data_ptr : np.float64*
        A pointer to the data of the numpy array.
    wscale : double
        The scale of the vector.
    n_features : int
        The number of features (= dimensionality of ``w``).
    sq_norm : double
        The squared norm of ``w``.
    """

    def __cinit__(self, np.ndarray[DOUBLE, ndim=1, mode='c'] w):
        if w.shape[0] > INT_MAX:
            raise ValueError("more than %d features not supported" % INT_MAX)
        self.w = w
        self.w_data_ptr = <DOUBLE *>w.data
        self.wscale = 1.0
        self.n_features = w.shape[0]
        self.sq_norm = np.dot(w, w)

    cdef void add(self, DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                  int xnnz, double c) nogil:
        """Scales example x by constant c and adds it to the weight vector.

        This operation updates ``sq_norm``.

        Parameters
        ----------
        x_data_ptr : double*
            The array which holds the feature values of ``x``.
        x_ind_ptr : np.int32*
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
        cdef DOUBLE* w_data_ptr = self.w_data_ptr

        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            innerprod += (w_data_ptr[idx] * val)
            xsqnorm += (val * val)
            w_data_ptr[idx] += val * (c / wscale)

        self.sq_norm += (xsqnorm * c * c) + (2.0 * innerprod * wscale * c)

    cdef double dot(self, DOUBLE *x_data_ptr, INTEGER *x_ind_ptr, int xnnz) nogil:
        """Computes the dot product of a sample x and the weight vector.

        Parameters
        ----------
        x_data_ptr : double*
            The array which holds the feature values of ``x``.
        x_ind_ptr : np.int32*
            The array which holds the feature indices of ``x``.
        xnnz : int
            The number of non-zero features of ``x``.

        Returns
        -------
        innerprod : double
            The inner product of ``x`` and ``w``.
        """
        cdef int j
        cdef int idx
        cdef double innerprod = 0.0
        cdef DOUBLE* w_data_ptr = self.w_data_ptr
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
            with gil:
                self.reset_wscale()

    cdef void reset_wscale(self) nogil:
        """Scales each coef of ``w`` by ``wscale`` and resets it to 1. """
        # self.w *= self.wscale, but without the GIL.
        dscal(<int>self.w.shape[0], self.wscale, <double *>self.w.data, 1)
        self.wscale = 1.0

    cdef double norm(self) nogil:
        """The L2 norm of the weight vector. """
        return sqrt(self.sq_norm)
