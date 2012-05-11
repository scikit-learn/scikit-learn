# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.


import numpy as np


cimport numpy as np
cimport cython


cdef class WeightVector(object):
    """Parameter vector for a generalized linear model.

    The model parameters consits of a weight vector and an
    intercept term (= a constant offset).
    The weight vector is represented by a scalar and a ndarray.
    The intercept term is a ndarray.

    The class provides methods to ``add`` a sparse vector
    and scale the vector.
    Representing a vector explicitly as a scalar times a
    vector allows for efficient scaling operations.

    Attributes
    ----------
    w : ndarray, dtype=np.float64, order='C', shape=(n_features, K)
        The numpy array which backs the weight vector.
    w_data_ptr : np.float64*
        A pointer to the data of the ndarray.
    w_scale : double
        The scale of the vector.
    intercept : ndarray, dtype=np.float64, order='C', shape=(K,)
        The intercept terms.
    intercept_data_ptr : np.float64*
        A pointer to the data of the ndarray
    n_features : int
        The number of features (= dimensionality of ``w``).
    sq_norm : double
        The squared norm of ``w``.
    K : Pysize_t
        The number of classes; 1 for regression and binary classification
        ``n_classes`` for multi-class classification.
    fit_intercept : int
        Whether or not to update the intercept term
    intercept_decay : double
        Scaling factor for the update of the intercept term.
    """

    def __cinit__(self, np.ndarray[np.float64_t, ndim=2] w,
                  np.ndarray[np.float64_t, ndim=1, mode='c'] intercept,
                  int fit_intercept, double intercept_decay):
        self.w = w
        self.w_data_ptr = <np.float64_t *>w.data
        self.w_scale = 1.0
        self.n_features = w.shape[0]
        self.K = w.shape[1]
        self.sq_norm = np.dot(w.ravel(), w.ravel())
        self.intercept = intercept
        self.intercept_data_ptr = <np.float64_t *>intercept.data
        self.fit_intercept = fit_intercept
        self.intercept_decay = intercept_decay

    cdef void add(self, DTYPE *x_data_ptr, INTEGER *x_ind_ptr,
                  int xnnz, int k, double c):
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
        k : int
            The class indices to compute the offset in the vector.
        c : double
            The scaling constant for the example.
        """
        cdef int j
        cdef int idx
        cdef Py_ssize_t K = self.K
        cdef double val
        cdef double innerprod = 0.0
        cdef double xsqnorm = 0.0

        # the next two lines save a factor of 2!
        cdef double w_scale = self.w_scale
        cdef np.float64_t* w_data_ptr = self.w_data_ptr

        for j in range(xnnz):
            if K > 1:
                idx = (x_ind_ptr[j] * K + k)
            else:
                idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            innerprod += (w_data_ptr[idx] * val)
            xsqnorm += (val * val)
            w_data_ptr[idx] += val * (c / w_scale)

        if self.fit_intercept:
            self.intercept_data_ptr[k] += c * self.intercept_decay

        self.sq_norm += (xsqnorm * c * c) + (2.0 * innerprod * w_scale * c)

    cdef void dot(self, DTYPE *x_data_ptr, INTEGER *x_ind_ptr, int xnnz, double *out):
        """Computes the dot product of a sample x and the weight vector.

        Parameters
        ----------
        x_data_ptr : double*
            The array which holds the feature values of ``x``.
        x_ind_ptr : np.int32*
            The array which holds the feature indices of ``x``.
        xnnz : int
            The number of non-zero features of ``x``.
        out : double*
            An array of length ``K`` in which the results are stored.
        """
        cdef int j
        cdef int idx
        cdef int k
        cdef DTYPE val
        cdef int class_offset = 0
        cdef np.float64_t* w_data_ptr # = self.w_data_ptr
        cdef Py_ssize_t K = self.K
        cdef Py_ssize_t n_features = self.n_features
        cdef double w_scale = self.w_scale

        for k in range(K):
            out[k] = 0.0

        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            w_data_ptr = self.w_data_ptr + (idx * K)
            for k in range(K):
                out[k] += w_data_ptr[k] * val

        for k in range(K):
            out[k] *= w_scale
            out[k] += self.intercept_data_ptr[k]

    cdef void scale(self, double c):
        """Scales the weight vector by a constant ``c``.

        It updates ``w_scale`` and ``sq_norm``. If ``w_scale`` gets too
        small we call ``reset_scale``."""
        self.w_scale *= c
        self.sq_norm *= (c * c)
        if self.w_scale < 1e-9:
            self.reset_scale()

    cdef void reset_scale(self):
        """Scales each coef of ``w`` by ``w_scale`` and resets it to 1. """
        self.w *= self.w_scale
        self.w_scale = 1.0

    cdef double norm(self):
        """The L2 norm of the weight vector. """
        return sqrt(self.sq_norm)

    def to_array(self):
        if self.w_scale != 1.0:
            self.reset_scale()
        return self.w, self.intercept
