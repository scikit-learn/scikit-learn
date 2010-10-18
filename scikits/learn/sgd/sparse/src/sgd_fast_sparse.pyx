# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
from __future__ import division

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
ctypedef np.float32_t FLOAT
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
    cpdef double loss(self,double p,double y):
        raise NotImplementedError()
    cpdef double dloss(self,double p,double y):
        raise NotImplementedError()


cdef class Classification(LossFunction):
    """Base class for loss functions for classification."""
    cpdef double loss(self,double p,double y):
        raise NotImplementedError()
    cpdef double dloss(self,double p,double y):
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
    cpdef double loss(self,double p,double y):
        cdef double z = p*y
        if z >= 1:
            return 0
        elif z >= -1:
            return (1-z) * (1-z) 
        else:
            return -4*z

    cpdef double dloss(self,double p,double y):
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
    cpdef  double loss(self,double p,double y):
        cdef double z = p*y
        if z < 1.0:
            return (1 - z)
        return 0
    cpdef  double dloss(self,double p,double y):
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
    cpdef double loss(self,double p,double y):
        cdef double z = p*y
        if z > 18:
            return exp(-z)
        if z < -18:
            return -z * y
        return log(1.0+exp(-z)) 

    cpdef  double dloss(self,double p,double y):
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
    cpdef  double loss(self,double p,double y):
        return 0.5 * (p-y) * (p-y)
    cpdef  double dloss(self,double p,double y):
        return y - p

    def __reduce__(self):
        return SquaredError,()

cdef class Huber(Regression):
    """
    """
    cdef double c
    def __init__(self,c):
        self.c = c
    cpdef  double loss(self,double p,double y):
        cdef double r = p-y
        cdef double abs_r = abs(r)
        if abs_r <= self.c:
            return 0.5 * r * r
        else:
            return self.c * abs_r - (0.5*self.c*self.c)

    cpdef  double dloss(self,double p,double y):
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
def plain_sgd(np.ndarray[DOUBLE, ndim=1] w,
              DOUBLE intercept,
              LossFunction loss,
              int penalty_type,
              double alpha, double rho,
              np.ndarray[FLOAT, ndim=1] X_data,
              np.ndarray[INTEGER, ndim=1] X_indices,
              np.ndarray[INTEGER, ndim=1] X_indptr,
              np.ndarray[DOUBLE, ndim=1] Y,
              int n_iter, int fit_intercept,
              int verbose, int shuffle):
    """Cython implementation of SGD with different loss functions and
    penalties.
    
    """
    # get the data information into easy vars
    cdef unsigned int n_samples = Y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    # FIXME double to DOUBLE
    cdef DOUBLE *w_data_ptr = <DOUBLE *>w.data
    cdef FLOAT *X_data_ptr = <FLOAT *>X_data.data
    cdef INTEGER *X_indptr_ptr = <INTEGER *>X_indptr.data
    cdef INTEGER *X_indices_ptr = <INTEGER *>X_indices.data

    # FIXME unsined int?
    cdef np.ndarray[INTEGER, ndim=1, mode="c"] index = np.arange(n_samples,
                                                                 dtype = np.int32)
    cdef INTEGER *index_ptr = <INTEGER *>index.data
    cdef INTEGER offset = 0
    cdef INTEGER xnnz = 0
    cdef DOUBLE wscale = 1.0
    cdef DOUBLE eta = 0.0
    cdef DOUBLE p = 0.0
    cdef DOUBLE update = 0.0
    cdef DOUBLE sumloss = 0.0
    cdef DOUBLE wnorm = 0.0
    cdef DOUBLE t = 0.0
    cdef DOUBLE y = 0.0
    cdef unsigned int count = 0
    cdef unsigned int epoch = 0
    cdef unsigned int i = 0
    cdef unsigned int sample_idx = 0
    cdef np.ndarray[DOUBLE, ndim=1, mode="c"] q = None
    cdef DOUBLE *q_ptr
    if penalty_type != L2:
        q = np.zeros((n_features,), dtype = np.float64, order = "c")
        q_ptr = <DOUBLE *> q.data
    cdef double u = 0.0
    # computing eta0
    cdef DOUBLE typw = sqrt(1.0 / sqrt(alpha))
    cdef DOUBLE eta0 = typw /max(1.0, loss.dloss(-typw, 1.0))
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
                wscale *= (1 - rho * eta * alpha)
                if wscale < 1e-9:
                    w *= wscale
                    wscale = 1.0
            if penalty_type == L2 or penalty_type == ELASTICNET:
                u += ((1.0 - rho) * eta * alpha)
                #l1penalty(wscale, wdata, qdata, xdata, xnnz, u)
            t += 1
            count += 1
        # report epoche information
        if verbose > 0:
            wnorm = sqrt(np.dot(w, w) * wscale * wscale)
            print("Norm: %.2f, NNZs: %d, Bias: %.6f, T: %d, Avg. loss: %.6f" % (wnorm, w.nonzero()[0].shape[0], intercept, count, sumloss / count))
            print("Total training time: %.2f seconds." % (time()-t_start))

        # floating-point under-/overflow check.
        if np.any(np.isinf(w)) or np.any(np.isnan(w)) or np.isnan(intercept) or np.isinf(intercept):
            raise ValueError("floating-point under-/overflow occured.")     

    w *= wscale
    return w, intercept
    

## # ----------------------------------------
## # Python function for external prediction
## # ----------------------------------------
## def predict(np.ndarray x, np.ndarray w,
##             double bias):
##     """Computes x*w + b efficiently.

##     :arg x: the instance represented as a sparse vector. 
##     :type x: np.ndarray(dtype=bolt.sparsedtype)
##     :arg w: the weight vector represented as a dense vector.
##     :type w: np.ndarray(dtype=bolt.densedtype)
##     :arg b: the bias term (aka offset or intercept).
##     :type b: float
##     :returns: A double representing `x*w + b`.
##     """
##     cdef int xnnz = x.shape[0]
##     cdef int wdim = w.shape[0]
##     cdef double y = 0.0
##     if xnnz == 0:
##         y = bias
##     else:
##         y = dot_checked(<double *>w.data,<Pair *>x.data,xnnz,wdim) + bias
##     return y

cdef inline double max(DOUBLE a, DOUBLE b):
    return a if a >= b else b

cdef inline DOUBLE min(DOUBLE a, DOUBLE b):
    return a if a <= b else b

cdef DOUBLE dot(DOUBLE *w_data_ptr, FLOAT *X_data_ptr, INTEGER *X_indices_ptr,
                INTEGER offset, INTEGER xnnz):
    cdef DOUBLE sum = 0.0
    cdef int j
    for j from 0 <= j < xnnz:
        sum += w_data_ptr[X_indices_ptr[offset + j]] * X_data_ptr[offset + j]
    return sum

cdef DOUBLE add(DOUBLE *w_data_ptr, DOUBLE wscale, FLOAT *X_data_ptr,
                INTEGER *X_indices_ptr, INTEGER offset, INTEGER xnnz, DOUBLE c):
    """Scales example x by constant c and adds it to the weight vector w. 
    """
    cdef INTEGER j
    cdef INTEGER idx
    cdef DOUBLE val
    cdef DOUBLE innerprod = 0.0
    cdef DOUBLE xsqnorm = 0.0
    for j from 0 <= j < xnnz:
        idx = X_indices_ptr[offset + j]
        val = X_data_ptr[offset + j]
        innerprod += (w_data_ptr[idx] * val)
        xsqnorm += (val * val)
        w_data_ptr[idx] += val * (c / wscale)    
    return (xsqnorm * c * c) + (2.0 * innerprod * wscale * c)

## cdef double dot_checked(double *w, Pair *x, int nnz, int wdim):
##     """ Checked version of dot product. Ignores features in x
##     with a higher index than dimension of w. 
##     """
##     cdef double sum = 0.0
##     cdef Pair pair
##     cdef int i
##     for i from 0 <= i < nnz:
##         pair = x[i]
##         if pair.idx < wdim:
##             sum +=w[pair.idx]*pair.val
##     return sum

## cdef double add(double *w, double wscale, Pair *x, int nnz, double c):
##     """Scales example x by constant c and adds it to the weight vector w. 
##     """
##     cdef Pair pair
##     cdef int i
##     cdef double innerprod = 0.0
##     cdef double xsqnorm = 0.0
##     for i from 0 <= i < nnz:
##         pair = x[i]
##         innerprod += (w[pair.idx] * pair.val)
##         xsqnorm += (pair.val*pair.val)
##         w[pair.idx] += pair.val * (c / wscale)
        
##     return (xsqnorm * c * c) + (2.0 * innerprod * wscale * c)

## # ----------------------------------------
## # Extension type for Stochastic Gradient Descent
## # ----------------------------------------

## cdef class SGD:

##     cdef int epochs
##     cdef double reg
##     cdef LossFunction loss
##     cdef int norm
##     cdef double alpha
    
##     def __init__(self, loss, reg, epochs = 5, norm = 2, alpha = 0.85):
##         """

##         :arg loss: The :class:`LossFunction` (default ModifiedHuber) .
##         :arg reg: The regularization parameter lambda (>0).
##         :type reg: float.
##         :arg epochs: The number of iterations through the dataset.
##         :type epochs: int
##         :arg norm: Whether to minimize the L1, L2 norm or the Elastic Net.
##         :type norm: 1 or 2 or 3
##         :arg alpha: The elastic net penality parameter. A value of 1 amounts to L2 regularization whereas a value of 0 gives L1 penalty. 
##         :type alpha: float (0 <= alpha <= 1)
##         """
##         if loss == None:
##             raise ValueError("Loss function must not be None.")
##         if reg < 0.0:
##             raise ValueError("reg must be larger than 0. ")
##         if norm not in [1,2,3]:
##             raise ValueError("norm must be in {1,2,3}. ")
##         if alpha > 1.0 or alpha < 0.0:
##             raise ValueError("alpha must be in [0,1]. ")
##         self.loss = loss
##         self.reg = reg
##         self.epochs = epochs
##         self.norm = norm
##         self.alpha = alpha

##     def __reduce__(self):
##         return SGD,(self.loss,self.reg, self.epochs, self.norm, self.alpha)

##     def train(self, model, dataset, verbose = 0, shuffle = False):
##         """Train `model` on the `dataset` using SGD.

##         :arg model: The :class:`bolt.model.LinearModel` that is going to be trained. 
##         :arg dataset: The :class:`bolt.io.Dataset`. 
##         :arg verbose: The verbosity level. If 0 no output to stdout.
##         :arg shuffle: Whether or not the training data should be shuffled after each epoch. 
##         """
##         self._train(model, dataset, verbose, shuffle)

##     cdef void _train(self,model, dataset, verbose, shuffle):
        
##         cdef LossFunction loss = self.loss
##         cdef int m = model.m
##         cdef int n = dataset.n
##         cdef double reg = self.reg

##         cdef np.ndarray[np.float64_t, ndim=1, mode="c"] w = model.w
##         # weight vector w as c array
##         cdef double *wdata = <double *>w.data
##         # the scale of w
##         cdef double wscale = 1.0

##         # training instance
##         cdef np.ndarray x = None
##         cdef Pair *xdata = NULL

##         cdef double y = 0.0
        
##         # Variables for penalty term
##         cdef int norm = self.norm
##         cdef np.ndarray[np.float64_t, ndim=1, mode="c"] q = None
##         cdef double *qdata = NULL
##         cdef double u = 0.0
##         if norm == 1 or norm == 3:
##             q = np.zeros((m,), dtype = np.float64, order = "c" )
##             qdata = <double *>q.data
            
##         cdef double alpha = 1.0
##         if norm == 1:
##             alpha = 0.0
##         elif norm == 3:
##             alpha = self.alpha

##         # bias term (aka offset or intercept)
##         cdef int usebias = 1
##         if model.biasterm == False:
##             usebias = 0

##         cdef double b = 0.0,p = 0.0, wnorm = 0.0, t = 0.0, update = 0.0,sumloss = 0.0, eta = 0.0
##         cdef int xnnz = 0, count = 0, i = 0, e = 0
        
##         # computing eta0
##         cdef double typw = sqrt(1.0 / sqrt(reg))
##         cdef double eta0 = typw /max(1.0,loss.dloss(-typw,1.0))
##         t = 1.0 / (eta0 * reg)
##         t1=time()
##         for e from 0 <= e < self.epochs:
##             if verbose > 0:
##                 print("-- Epoch %d" % (e+1))
##             if shuffle:
##                 dataset.shuffle()
##             for x,y in dataset:
##                 eta = 1.0 / (reg * t)
##                 xnnz = x.shape[0]
##                 xdata = <Pair *>x.data
##                 p = (dot(wdata, xdata, xnnz) * wscale) + b
##                 sumloss += loss.loss(p,y)
##                 update = eta * loss.dloss(p,y)
##                 if update != 0:
##                     add(wdata, wscale, xdata,
##                         xnnz,update)
##                     if usebias == 1:
##                         b += update * 0.01

##                 if norm != 1:
##                     wscale *= (1 - alpha * eta * reg)
##                     if wscale < 1e-9:
##                         w*=wscale
##                         wscale = 1
##                 if norm == 1 or norm == 3:
##                     u += ((1-alpha) * eta * reg)
##                     l1penalty(wscale, wdata, qdata, xdata, xnnz, u)
                
##                 t += 1
##                 count += 1

##             # report epoche information
##             if verbose > 0:
##                 wnorm = sqrt(np.dot(w,w) * wscale * wscale)
##                 print("Norm: %.2f, NNZs: %d, Bias: %.6f, T: %d, Avg. loss: %.6f" % (wnorm,w.nonzero()[0].shape[0],b,count,sumloss/count))
##                 print("Total training time: %.2f seconds." % (time()-t1))

##         # floating-point under-/overflow check.
##         if np.any(np.isinf(w)) or np.any(np.isnan(w))or np.isnan(b) or np.isinf(b):
##             raise ValueError("floating-point under-/overflow occured.")
##         if norm == 3:
##             # FIXME rescale naive elastic net coefficient?
##             model.w = w * wscale #* (1.0 + alpha)
##         else:
##             model.w = w * wscale
##         model.bias = b

## cdef void l1penalty(double wscale, double *w, double *q,
##                     Pair *x, int nnz, double u):
##     cdef double z = 0.0
##     cdef Pair pair
##     cdef int i,j
##     for i from 0 <= i < nnz:
##         pair = x[i]
##         j = pair.idx
##         z = w[j]
##         if (wscale * w[j]) > 0:
##             w[j] = max(0,w[j] - ((u + q[j])/wscale) )
##         elif (wscale * w[j]) < 0:
##             w[j] = min(0,w[j] + ((u - q[j])/wscale) )
##         q[j] += (wscale * (w[j] - z))

