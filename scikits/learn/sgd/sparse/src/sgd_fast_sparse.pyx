# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.

import numpy as np
import sys
from time import time

cimport numpy as np
cimport cython

cdef extern from "math.h":
    cdef extern double exp(double x)
    cdef extern double log(double x)
    cdef extern double sqrt(double x)

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER

DEF L1 = 1
DEF L2 = 2
DEF ELASTICNET = 3

# ----------------------------------------
# Extension Types for Loss Functions
# ----------------------------------------

cdef class LossFunction:
    """Base class for convex loss functions"""
    cpdef double loss(self, double p, double y):
        """Evaluate the loss function.
        
        :arg p: The prediction.
        :type p: double
        :arg y: The true value.
        :type y: double
        :returns: double"""
        raise NotImplementedError()
    cpdef double dloss(self, double p, double y):
        """Evaluate the derivative of the loss function.
        
        :arg p: The prediction.
        :type p: double
        :arg y: The true value.
        :type y: double
        :returns: double"""
        raise NotImplementedError()

cdef class Regression(LossFunction):
    """Base class for loss functions for regression."""
    cpdef double loss(self,double p, double y):
        raise NotImplementedError()
    cpdef double dloss(self,double p, double y):
        raise NotImplementedError()


cdef class Classification(LossFunction):
    """Base class for loss functions for classification."""
    cpdef double loss(self, double p, double y):
        raise NotImplementedError()
    cpdef double dloss(self, double p, double y):
        raise NotImplementedError()

cdef class ModifiedHuber(Classification):
    """Modified Huber loss function for binary
    classification tasks with y in {-1,1}.
    Its equivalent to quadratically smoothed SVM
    with gamma = 2.

    See T. Zhang 'Solving
    Large Scale Linear Prediction Problems Using
    Stochastic Gradient Descent', ICML'04.
    """
    cpdef double loss(self, double p,double y):
        cdef double z = p*y
        if z >= 1:
            return 0
        elif z >= -1:
            return (1-z) * (1-z) 
        else:
            return -4*z

    cpdef double dloss(self, double p, double y):
        cdef double z = p*y
        if z >= 1:
            return 0
        elif z >= -1:
            return 2*(1-z)*y
        else:
            return 4*y

    def __reduce__(self):
        return ModifiedHuber,()

cdef class Hinge(Classification):
    """SVM classification loss for binary
    classification tasks with y in {-1,1}.
    """
    cpdef  double loss(self, double p, double y):
        cdef double z = p*y
        if z < 1.0:
            return (1 - z)
        return 0
    cpdef  double dloss(self, double p, double y):
        cdef double z = p*y
        if z < 1.0:
            return y
        return 0

    def __reduce__(self):
        return Hinge,()


cdef class Log(Classification):
    """Logistic regression loss for binary classification
    tasks with y in {-1,1}.
    """
    cpdef double loss(self, double p, double y):
        cdef double z = p*y
        if z > 18:
            return exp(-z)
        if z < -18:
            return -z * y
        return log(1.0+exp(-z)) 

    cpdef  double dloss(self, double p, double y):
        cdef double z = p*y
        if z > 18:
            return exp(-z) * y
        if z < -18:
            return y
        return y / (exp(z) + 1.0)

    def __reduce__(self):
        return Log,()

cdef class SquaredError(Regression):
    """
    """
    cpdef  double loss(self, double p, double y):
        return 0.5 * (p-y) * (p-y)
    cpdef  double dloss(self, double p, double y):
        return y - p

    def __reduce__(self):
        return SquaredError,()

cdef class Huber(Regression):
    """
    """
    cdef double c
    def __init__(self,c):
        self.c = c
    cpdef  double loss(self, double p, double y):
        cdef double r = p-y
        cdef double abs_r = abs(r)
        if abs_r <= self.c:
            return 0.5 * r * r
        else:
            return self.c * abs_r - (0.5*self.c*self.c)

    cpdef  double dloss(self, double p, double y):
        cdef double r = y - p 
        cdef double abs_r = abs(r)
        if abs_r <= self.c:
            return r
        elif r > 0:
            return self.c
        else:
            return -self.c

    def __reduce__(self):
        return Huber,(self.c,)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def plain_sgd(np.ndarray[double, ndim=1] w,
              double intercept,
              LossFunction loss,
              int penalty_type,
              double alpha, double rho,
              np.ndarray[double, ndim=1] X_data,
              np.ndarray[int, ndim=1] X_indices,
              np.ndarray[int, ndim=1] X_indptr,
              np.ndarray[double, ndim=1] Y,
              int n_iter, int fit_intercept,
              int verbose, int shuffle):
    """Cython implementation of SGD with different loss functions and
    penalties.
    
    """
    # get the data information into easy vars
    cdef unsigned int n_samples = Y.shape[0]
    cdef unsigned int n_features = w.shape[0]
    
    cdef double *w_data_ptr = <double *>w.data
    cdef double *X_data_ptr = <double *>X_data.data
    cdef int *X_indptr_ptr = <int *>X_indptr.data
    cdef int *X_indices_ptr = <int *>X_indices.data

    # FIXME unsined int?
    cdef np.ndarray[int, ndim=1, mode="c"] index = np.arange(n_samples,
                                                             dtype = np.int32)
    cdef int *index_ptr = <int *>index.data
    cdef int offset = 0
    cdef int xnnz = 0
    cdef double wscale = 1.0
    cdef double eta = 0.0
    cdef double p = 0.0
    cdef double update = 0.0
    cdef double sumloss = 0.0
    cdef double wnorm = 0.0
    cdef double t = 0.0
    cdef double y = 0.0
    cdef unsigned int count = 0
    cdef unsigned int epoch = 0
    cdef unsigned int i = 0
    cdef unsigned int sample_idx = 0
    cdef np.ndarray[double, ndim=1, mode="c"] q = None
    cdef double *q_data_ptr
    if penalty_type != L2:
        q = np.zeros((n_features,), dtype = np.float64, order = "c")
        q_data_ptr = <double *> q.data
    cdef double u = 0.0
    # computing eta0
    cdef double typw = sqrt(1.0 / sqrt(alpha))
    cdef double eta0 = typw / max(1.0, loss.dloss(-typw, 1.0))
    t = 1.0 / (eta0 * alpha)
    t_start = time()
    for epoch from 0 <= epoch < n_iter:
        if verbose > 0:
            print("-- Epoch %d" % (epoch + 1))
        if shuffle:
            np.random.shuffle(index)
        for i from 0 <= i < n_samples:
            sample_idx = index[i]
            offset = X_indptr_ptr[sample_idx]
            xnnz = X_indptr_ptr[sample_idx + 1] - offset
            y = Y[sample_idx]
            eta = 1.0 / (alpha * t)
            p = (dot(w_data_ptr, X_data_ptr, X_indices_ptr,
                     offset, xnnz) * wscale) + intercept
            sumloss += loss.loss(p, y)
            update = eta * loss.dloss(p, y)
            if update != 0:
                add(w_data_ptr, wscale, X_data_ptr, X_indices_ptr,
                    offset, xnnz, update)
                if fit_intercept == 1:
                    intercept += update * 0.01
            if penalty_type != L1:
                wscale *= (1 - (rho * eta * alpha))
                if wscale < 1e-9:
                    w *= wscale
                    wscale = 1.0
            if penalty_type == L1 or penalty_type == ELASTICNET:
                u += ((1 - rho) * eta * alpha)
                l1penalty(w_data_ptr, wscale, q_data_ptr,
                          X_indices_ptr, offset, xnnz, u)
            t += 1
            count += 1
        if penalty_type == L1 or penalty_type == ELASTICNET:
            u += ((1 - rho) * eta * alpha)
            finall1penalty(w_data_ptr, wscale, n_features, q_data_ptr, u)
	    
        # report epoche information
        if verbose > 0:
            wnorm = sqrt(np.dot(w, w) * wscale * wscale)
            print("Norm: %.2f, NNZs: %d, "\
            "Bias: %.6f, T: %d, Avg. loss: %.6f" % (wnorm, 
                                                    w.nonzero()[0].shape[0],
                                                    intercept, count, 
                                                    sumloss / count))
            print("Total training time: %.2f seconds." % (time()-t_start))

        # floating-point under-/overflow check.
        if np.any(np.isinf(w)) or np.any(np.isnan(w)) \
           or np.isnan(intercept) or np.isinf(intercept):
            raise ValueError("floating-point under-/overflow occured.")     

    w *= wscale
    return w, intercept


cdef inline double max(double a, double b):
    return a if a >= b else b

cdef inline double min(double a, double b):
    return a if a <= b else b

cdef double dot(double *w_data_ptr, double *X_data_ptr, int *X_indices_ptr,
                int offset, int xnnz):
    cdef double sum = 0.0
    cdef int j
    for j from 0 <= j < xnnz:
        sum += w_data_ptr[X_indices_ptr[offset + j]] * X_data_ptr[offset + j]
    return sum

cdef double add(double *w_data_ptr, double wscale, double *X_data_ptr,
                int *X_indices_ptr, int offset, int xnnz, double c):
    """Scales example x by constant c and adds it to the weight vector w. 
    """
    cdef int j
    cdef int idx
    cdef double val
    cdef double innerprod = 0.0
    cdef double xsqnorm = 0.0
    for j from 0 <= j < xnnz:
        idx = X_indices_ptr[offset + j]
        val = X_data_ptr[offset + j]
        innerprod += (w_data_ptr[idx] * val)
        xsqnorm += (val * val)
        w_data_ptr[idx] += val * (c / wscale)    
    return (xsqnorm * c * c) + (2.0 * innerprod * wscale * c)

cdef void l1penalty(double *w_data_ptr, double wscale, double *q_data_ptr,
                    int *X_indices_ptr, int offset, int xnnz, double u):
    """Applys the L1 penalty to each updated feature. 
    This implements the truncated gradient approach by 
    [Tsuruoka, Y., Tsujii, J., and Ananiadou, S., 2009].
    """
    cdef double z = 0.0
    cdef int j = 0
    cdef int idx = 0
    for j from 0 <= j < xnnz:
        idx = X_indices_ptr[offset + j]
        z = w_data_ptr[idx]
        if (wscale * w_data_ptr[idx]) > 0:
            w_data_ptr[idx] = max(0, w_data_ptr[idx] - ((u + q_data_ptr[idx])
                                                        / wscale) )
        elif (wscale * w_data_ptr[idx]) < 0:
            w_data_ptr[idx] = min(0, w_data_ptr[idx] + ((u - q_data_ptr[idx])
                                                        / wscale) )
        q_data_ptr[idx] += (wscale * (w_data_ptr[idx] - z))

cdef void finall1penalty(double *w_data_ptr, double wscale, 
                         unsigned int n_features, 
                         double *q_data_ptr, double u):
    """Applys the L1 penalty to all feature. 
    This implements the truncated gradient approach by 
    [Tsuruoka, Y., Tsujii, J., and Ananiadou, S., 2009].

    Experimental: this was proposed by Bob Carpenter (LingPipe). 
    
    """
    cdef double z = 0.0
    cdef int j = 0
    for j from 0 <= j < n_features:
        z = w_data_ptr[j]
        if (wscale * w_data_ptr[j]) > 0:
            w_data_ptr[j] = max(0, w_data_ptr[j] - ((u + q_data_ptr[j])
                                                    / wscale) )
        elif (wscale * w_data_ptr[j]) < 0:
            w_data_ptr[j] = min(0, w_data_ptr[j] + ((u - q_data_ptr[j])
                                                    / wscale) )
        q_data_ptr[j] += (wscale * (w_data_ptr[j] - z))
    
