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

    def __cinit__(self, np.ndarray[np.float64_t, ndim=2, mode='c'] w,
                  np.ndarray[np.float64_t, ndim=1, mode='c'] intercept,
                  int fit_intercept, double intercept_decay):
        self.w = w
        self.w_data_ptr = <np.float64_t *>w.data
        self.w_scale = 1.0
        self.n_features = w.shape[1]
        self.K = w.shape[0]
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
        cdef double val
        cdef double innerprod = 0.0
        cdef double xsqnorm = 0.0
        cdef int class_offset = 0
        if k > 0:
            class_offset = k * self.n_features

        # the next two lines save a factor of 2!
        cdef double w_scale = self.w_scale
        cdef np.float64_t* w_data_ptr = self.w_data_ptr + class_offset

        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            innerprod += (w_data_ptr[idx] * val)
            xsqnorm += (val * val)
            w_data_ptr[idx] += val * (c / w_scale)

        if self.fit_intercept:
            self.intercept_data_ptr[k] += c * self.intercept_decay

        self.sq_norm += (xsqnorm * c * c) + (2.0 * innerprod * w_scale * c)

    cdef double dot(self, DTYPE *x_data_ptr, INTEGER *x_ind_ptr, int xnnz, int k):
        """Computes the dot product of a sample x and the weight vector.

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

        Returns
        -------
        innerprod : double
            The inner product of ``x`` and ``w``.
        """
        cdef int j
        cdef int idx
        cdef double innerprod = 0.0
        cdef int class_offset = 0
        if k > 0:
            class_offset = k * self.n_features
        cdef np.float64_t* w_data_ptr = self.w_data_ptr + class_offset
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            innerprod += w_data_ptr[idx] * x_data_ptr[j]
        innerprod *= self.w_scale

        # add intercept
        innerprod += self.intercept_data_ptr[k]
        return innerprod

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



cdef class AveragedWeightVector(WeightVector):
    """A weight vector that is averaged over its mutations.

    The intercept is not subject to the averaging.

    Attributes
    ----------
    w : ndarray, dtype=np.float64, order='C', shape=(n_features, K)
        The numpy array which backs the weight vector.
    w_data_ptr : np.float64*
        A pointer to the data of the ndarray.
    w_scale : double
        The scale of the vector.
    w_bar : ndarray, dtype=np.float64, order='C', shape=(n_features, K)
        The averaged weight vector
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

    def __cinit__(self, np.ndarray[np.float64_t, ndim=2, mode='c'] w,
                  np.ndarray[np.float64_t, ndim=1, mode='c'] intercept,
                  int fit_intercept, double intercept_decay,
                  int n_updates):
        self.w_bar = w
        self.w = np.copy(w)
        self.w_data_ptr = <np.float64_t *>w.data
        self.w_bar_data_ptr = <np.float64_t *>self.w_bar.data
        self.w_scale = 1.0
        self.w_bar_scale = 1.0
        self.n_features = w.shape[1]
        self.K = w.shape[0]
        self.sq_norm = np.dot(w.ravel(), w.ravel())
        self.intercept = intercept
        self.intercept_data_ptr = <np.float64_t *>intercept.data
        self.fit_intercept = fit_intercept
        self.intercept_decay = intercept_decay
        self.n_updates = n_updates
        self.n_dots = 0

    cdef void add(self, DTYPE *x_data_ptr, INTEGER *x_ind_ptr,
                  int xnnz, int k, double c):
        """Scales example x by constant c and adds it to the weight vector.

        Each update is applied to both ``w`` and ``w_bar``.
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
        cdef double val
        cdef double innerprod = 0.0
        cdef double xsqnorm = 0.0
        cdef int class_offset = 0

        # the next two lines save a factor of 2!
        cdef double w_scale = self.w_scale
        cdef np.float64_t* w_data_ptr = self.w_data_ptr + class_offset
        cdef np.float64_t* w_bar_data_ptr = self.w_bar_data_ptr + class_offset

        # mutation will reside in w_bar for the next `u` updates
        cdef double u = (self.n_updates - (self.n_dots / self.K) + 1)

        if k > 0:
            class_offset = k * self.n_features

        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            innerprod += (w_data_ptr[idx] * val)
            xsqnorm += (val * val)
            w_data_ptr[idx] += val * c / w_scale
            w_bar_data_ptr[idx] += val * u * c / w_scale

        if self.fit_intercept:
            self.intercept_data_ptr[k] += c * self.intercept_decay

        self.sq_norm += (xsqnorm * c * c) + (2.0 * innerprod * w_scale * c)

    cdef double dot(self, DTYPE *x_data_ptr, INTEGER *x_ind_ptr, int xnnz, int k):
        """Computes the dot product of a sample x and the weight vector.

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

        Returns
        -------
        innerprod : double
            The inner product of ``x`` and ``w``.
        """
        cdef double innerprod = WeightVector.dot(self, x_data_ptr,
                                                 x_ind_ptr, xnnz, k)
        self.n_dots += 1
        return innerprod

    cdef void scale(self, double c):
        """Scales the weight vector by a constant ``c``.

        It updates ``w_scale`` and ``sq_norm``. If ``w_scale`` gets too
        small we call ``reset_scale``. """
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
        """Weight vector representation as ndarrays.

        Divides ``w_bar`` by the number of updates.
        """
        if self.w_scale != 1.0:
            self.reset_scale()
        w = self.w_bar / self.n_updates
        return w, self.intercept
