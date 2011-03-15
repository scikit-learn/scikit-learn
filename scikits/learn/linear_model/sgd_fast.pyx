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
    """Base class for loss functions for regression"""

    cpdef double loss(self,double p, double y):
        raise NotImplementedError()

    cpdef double dloss(self,double p, double y):
        raise NotImplementedError()


cdef class Classification(LossFunction):
    """Base class for loss functions for classification"""

    cpdef double loss(self, double p, double y):
        raise NotImplementedError()

    cpdef double dloss(self, double p, double y):
        raise NotImplementedError()


cdef class ModifiedHuber(Classification):
    """Modified Huber loss for binary classification with y in {-1, 1}

    This is equivalent to quadratically smoothed SVM with gamma = 2.

    See T. Zhang 'Solving Large Scale Linear Prediction Problems Using
    Stochastic Gradient Descent', ICML'04.
    """
    cpdef double loss(self, double p, double y):
        cdef double z = p * y
        if z >= 1.0:
            return 0.0
        elif z >= -1.0:
            return (1.0 - z) * (1.0 - z)
        else:
            return -4.0 * z

    cpdef double dloss(self, double p, double y):
        cdef double z = p * y
        if z >= 1.0:
            return 0.0
        elif z >= -1.0:
            return 2.0 * (1.0 - z) * y
        else:
            return 4.0 * y

    def __reduce__(self):
        return ModifiedHuber, ()


cdef class Hinge(Classification):
    """SVM loss for binary classification tasks with y in {-1,1}"""
    cpdef double loss(self, double p, double y):
        cdef double z = p * y
        if z < 1.0:
            return (1 - z)
        return 0.0

    cpdef double dloss(self, double p, double y):
        cdef double z = p * y
        if z < 1.0:
            return y
        return 0.0

    def __reduce__(self):
        return Hinge, ()


cdef class Log(Classification):
    """Logistic regression loss for binary classification with y in {-1, 1}"""

    cpdef double loss(self, double p, double y):
        cdef double z = p * y
        # approximately equal and saves the computation of the log
        if z > 18:
            return exp(-z)
        if z < -18:
            return -z * y
        return log(1.0 + exp(-z))

    cpdef double dloss(self, double p, double y):
        cdef double z = p * y
        # approximately equal and saves the computation of the log
        if z > 18.0:
            return exp(-z) * y
        if z < -18.0:
            return y
        return y / (exp(z) + 1.0)

    def __reduce__(self):
        return Log, ()


cdef class SquaredLoss(Regression):
    """Squared loss traditional used in linear regression."""
    cpdef double loss(self, double p, double y):
        return 0.5 * (p - y) * (p - y)

    cpdef double dloss(self, double p, double y):
        return y - p

    def __reduce__(self):
        return SquaredLoss, ()


cdef class Huber(Regression):
    """Huber regression loss

    Variant of the SquaredLoss that is robust to outliers (quadratic near zero,
    linear in for large errors).

    References
    ----------

    http://en.wikipedia.org/wiki/Huber_Loss_Function
    """

    def __init__(self,c):
        self.c = c

    cpdef double loss(self, double p, double y):
        cdef double r = p - y
        cdef double abs_r = abs(r)
        if abs_r <= self.c:
            return 0.5 * r * r
        else:
            return self.c * abs_r - (0.5 * self.c * self.c)

    cpdef double dloss(self, double p, double y):
        cdef double r = y - p
        cdef double abs_r = abs(r)
        if abs_r <= self.c:
            return r
        elif r > 0.0:
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
              np.ndarray[double, ndim=2] X,
              np.ndarray[double, ndim=1] Y,
              int n_iter, int fit_intercept,
              int verbose, int shuffle, int seed,
              double weight_pos, double weight_neg,
              np.ndarray[double, ndim=1] sample_weight):
    """Cython impl. of SGD for generic loss functions and penalties

    This implementation assumes X represented as a dense array of floats.

    Parameters
    ----------
    w : ndarray[double, ndim=1]
        The allocated coef_ vector.
    intercept : double
        The initial intercept
    loss : LossFunction
        A concrete LossFunction object.
    penalty_type : int
        The penalty 2 for L2, 1 for L1, and 3 for Elastic-Net.
    alpha : float
        The regularization parameter.
    rho : float
        The elastic net hyperparameter.
    X : ndarray[double, ndim=2]
        The dataset as a dense numpy array.
    Y : ndarray[double, ndim=1]
        The labels.
    n_iter : int
        The number of iterations (epochs).
    fit_intercept : int
        Whether or not to fit the intercept (1 or 0).
    verbose : int
        Print verbose output; 0 for quite.
    shuffle : int
        Whether to shuffle the training data before each epoch.
    weight_pos : float
        The weight of the positive class.
    weight_neg : float
        The weight of the negative class. 
    seed : int
        The seed of the pseudo random number generator to use when
        shuffling the data
    sample_weight : array, shape = [n_samples]
        The importance weight of each sample.

    Returns
    -------
    w : array, shape [n_features]
        The fitted weight vector. 
    intercept : float
        The fitted intercept term. 

    
    """

    # get the data information into easy vars
    cdef unsigned int n_samples = Y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    # Array strides to get to next feature or example
    cdef int row_stride = X.strides[0]
    cdef int elem_stride = X.strides[1]

    cdef double *w_data_ptr = <double *>w.data
    cdef double *X_data_ptr = <double *>X.data
    cdef double *Y_data_ptr = <double *>Y.data

    cdef double *sample_weight_data = <double *>sample_weight.data

    # Use index array for fast shuffling
    cdef np.ndarray[int, ndim=1, mode="c"] index = np.arange(n_samples,
                                                             dtype = np.int32)
    cdef int *index_data_ptr = <int *>index.data

    # helper variable
    cdef int offset = 0
    cdef double wscale = 1.0
    cdef double eta = 0.0
    cdef double p = 0.0
    cdef double update = 0.0
    cdef double sumloss = 0.0
    cdef double wnorm = 0.0
    cdef double t = 0.0
    cdef double y = 0.0
    cdef double class_weight = 1.0
    cdef unsigned int count = 0
    cdef unsigned int epoch = 0
    cdef unsigned int i = 0
    cdef int sample_idx = 0

    # q vector is only used for L1 regularization
    cdef np.ndarray[double, ndim=1, mode="c"] q = None
    cdef double *q_data_ptr
    if penalty_type != L2:
        q = np.zeros((n_features,), dtype = np.float64, order = "c")
        q_data_ptr = <double *> q.data
    cdef double u = 0.0

    # computing eta0, the initial learning rate
    cdef double typw = sqrt(1.0 / sqrt(alpha))
    cdef double eta0 = typw / max(1.0, loss.dloss(-typw, 1.0))

    # initialize the 1 / t learning rate schedule from eta0
    t = 1.0 / (eta0 * alpha)
    t_start = time()
    for epoch from 0 <= epoch < n_iter:
        if verbose > 0:
            print("-- Epoch %d" % (epoch + 1))
        if shuffle:
            np.random.RandomState(seed).shuffle(index)
        for i from 0 <= i < n_samples:
            sample_idx = index_data_ptr[i]
            offset = row_stride * sample_idx / elem_stride # row offset in elem
            y = Y_data_ptr[sample_idx]
            eta = 1.0 / (alpha * t)
            p = (dot(w_data_ptr, X_data_ptr, offset, n_features) * wscale
                ) + intercept
            sumloss += loss.loss(p, y)
            if y > 0:
                class_weight = weight_pos
            else:
                class_weight = weight_neg
            update = eta * loss.dloss(p, y) * class_weight * \
                sample_weight_data[sample_idx]
            if update != 0.0:
                add(w_data_ptr, wscale, X_data_ptr, offset, n_features, update)
                if fit_intercept == 1:
                    intercept += update
            if penalty_type != L1:
                wscale *= (1.0 - (rho * eta * alpha))
                if wscale < 1e-9:
                    w *= wscale
                    wscale = 1.0
            if penalty_type == L1 or penalty_type == ELASTICNET:
                u += ((1.0 - rho) * eta * alpha)
                l1penalty(w_data_ptr, wscale, q_data_ptr, n_features, u)
            t += 1
            count += 1

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


cdef double dot(double *w_data_ptr, double *X_data_ptr,
                int offset, unsigned int n_features):
    cdef double sum = 0.0
    cdef int j
    for j from 0 <= j < n_features:
        sum += w_data_ptr[j] * X_data_ptr[offset + j]
    return sum


cdef double add(double *w_data_ptr, double wscale, double *X_data_ptr,
                int offset, unsigned int n_features, double c):
    """Scales example x by constant c and adds it to the weight vector w"""
    cdef int j
    cdef int idx
    cdef double val
    cdef double innerprod = 0.0
    cdef double xsqnorm = 0.0
    for j from 0 <= j < n_features:
        val = X_data_ptr[offset + j]
        innerprod += (w_data_ptr[j] * val)
        xsqnorm += (val * val)
        w_data_ptr[j] += val * (c / wscale)

    # TODO this is needed for PEGASOS only
    return (xsqnorm * c * c) + (2.0 * innerprod * wscale * c)


cdef void l1penalty(double *w_data_ptr, double wscale, double *q_data_ptr,
                    unsigned int n_features, double u):
    """Apply the L1 penalty to each updated feature

    This implements the truncated gradient approach by
    [Tsuruoka, Y., Tsujii, J., and Ananiadou, S., 2009].

    FIXME: apply penalty over all features or only non-zero?
    Empirical results look better this way...
    """
    cdef double z = 0.0
    cdef int j = 0
    cdef int idx = 0
    for j from 0 <= j < n_features:
        z = w_data_ptr[j]
        if (wscale * w_data_ptr[j]) > 0.0:
            w_data_ptr[j] = max(0.0, w_data_ptr[j] - ((u + q_data_ptr[j])
                                                        / wscale) )
        elif (wscale * w_data_ptr[j]) < 0.0:
            w_data_ptr[j] = min(0.0, w_data_ptr[j] + ((u - q_data_ptr[j])
                                                        / wscale) )
        q_data_ptr[j] += (wscale * (w_data_ptr[j] - z))


