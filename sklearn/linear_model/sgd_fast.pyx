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

from libc.stdlib cimport free, malloc

from sklearn.utils.weight_vector cimport WeightVector
from sklearn.utils.seq_dataset cimport SequentialDataset


cdef extern from "math.h":
    cdef extern double exp(double x)
    cdef extern double log(double x)
    cdef extern double sqrt(double x)
    cdef extern double pow(double x, double y)

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER


# Penalty constans
DEF NO_PENALTY = 0
DEF L1 = 1
DEF L2 = 2
DEF ELASTICNET = 3

# Learning rate constants
DEF CONSTANT = 1
DEF OPTIMAL = 2
DEF INVSCALING = 3

# ----------------------------------------
# Extension Types for Loss Functions
# ----------------------------------------

cdef class LossFunction:
    """Base class for convex loss functions"""

    cpdef double loss(self, double p, double y):
        """Evaluate the loss function.

        Parameters
        ----------
        p : double
            The prediction, p = w^T x
        y : double
            The true value (aka target)

        Returns
        -------
        double
            The loss evaluated at `p` and `y`.
        """
        raise NotImplementedError()

    cpdef double dloss(self, double p, double y):
        """Evaluate the derivative of the loss function with respect to
        the prediction `p`.

        Parameters
        ----------
        p : double
            The prediction, p = w^T x
        y : double
            The true value (aka target)
        Returns
        -------
        double
            The derivative of the loss function w.r.t. `p`.
        """
        raise NotImplementedError()


cdef class Regression(LossFunction):
    """Base class for loss functions for regression"""

    cpdef double loss(self, double p, double y):
        raise NotImplementedError()

    cpdef double dloss(self, double p, double y):
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
            return 2.0 * (1.0 - z) * -y
        else:
            return -4.0 * y

    def __reduce__(self):
        return ModifiedHuber, ()


cdef class Hinge(Classification):
    """Hinge loss for binary classification tasks with y in {-1,1}

    Parameters
    ----------

    threshold : float > 0.0
        Margin threshold. When threshold=1.0, one gets the loss used by SVM.
        When threshold=0.0, one gets the loss used by the Perceptron.
    """

    cdef double threshold

    def __init__(self, double threshold=1.0):
        self.threshold = threshold

    cpdef double loss(self, double p, double y):
        cdef double z = p * y
        if z <= self.threshold:
            return (self.threshold - z)
        return 0.0

    cpdef double dloss(self, double p, double y):
        cdef double z = p * y
        if z <= self.threshold:
            return -y
        return 0.0

    def __reduce__(self):
        return Hinge, (self.threshold,)


cdef class Log(Classification):
    """Logistic regression loss for binary classification with y in {-1, 1}"""

    cpdef double loss(self, double p, double y):
        cdef double z = p * y
        # approximately equal and saves the computation of the log
        if z > 18:
            return exp(-z)
        if z < -18:
            return -z
        return log(1.0 + exp(-z))

    cpdef double dloss(self, double p, double y):
        cdef double z = p * y
        # approximately equal and saves the computation of the log
        if z > 18.0:
            return exp(-z) * -y
        if z < -18.0:
            return -y
        return -y / (exp(z) + 1.0)

    def __reduce__(self):
        return Log, ()


cdef class SquaredLoss(Regression):
    """Squared loss traditional used in linear regression."""
    cpdef double loss(self, double p, double y):
        return 0.5 * (p - y) * (p - y)

    cpdef double dloss(self, double p, double y):
        return p - y

    def __reduce__(self):
        return SquaredLoss, ()


cdef class Huber(Regression):
    """Huber regression loss

    Variant of the SquaredLoss that is robust to outliers (quadratic near zero,
    linear in for large errors).

    http://en.wikipedia.org/wiki/Huber_Loss_Function
    """

    cdef double c

    def __init__(self, double c):
        self.c = c

    cpdef double loss(self, double p, double y):
        cdef double r = p - y
        cdef double abs_r = abs(r)
        if abs_r <= self.c:
            return 0.5 * r * r
        else:
            return self.c * abs_r - (0.5 * self.c * self.c)

    cpdef double dloss(self, double p, double y):
        cdef double r = p - y
        cdef double abs_r = abs(r)
        if abs_r <= self.c:
            return r
        elif r > 0.0:
            return self.c
        else:
            return -self.c

    def __reduce__(self):
        return Huber, (self.c,)


def plain_sgd(np.ndarray[double, ndim=2, mode='c'] weights,
              np.ndarray[double, ndim=1, mode='c'] intercept,
              LossFunction loss,
              int penalty_type,
              double alpha, double rho,
              SequentialDataset dataset,
              int n_iter, int fit_intercept,
              int verbose, int shuffle, int seed,
              np.ndarray[double, ndim=1, mode='c'] class_weight,
              int learning_rate, double eta0,
              double power_t,
              double t=1.0,
              double intercept_decay=1.0):
    """Plain SGD for generic loss functions and penalties.

    Parameters
    ----------
    weights : ndarray[double, ndim=2], shape=(n_features, K)
        The allocated coef_ vector.
    intercept : ndarray[double, ndim=1], shape=(K,)
        The array of intercepts; shape
    loss : LossFunction
        A concrete ``LossFunction`` object.
    penalty_type : int
        The penalty 2 for L2, 1 for L1, and 3 for Elastic-Net.
    alpha : float
        The regularization parameter.
    rho : float
        The elastic net hyperparameter.
    dataset : SequentialDataset
        A concrete ``SequentialDataset`` object.
    n_iter : int
        The number of iterations (epochs).
    fit_intercept : int
        Whether or not to fit the intercept (1 or 0).
    verbose : int
        Print verbose output; 0 for quite.
    shuffle : int
        Whether to shuffle the training data before each epoch.
    class_weight : ndarray[double, ndim=1], shape=(n_classes,)
        The class weights.
    seed : int
        The seed of the pseudo random number generator to use when
        shuffling the data
    learning_rate : int
        The learning rate:
        (1) constant, eta = eta0
        (2) optimal, eta = 1.0/(t+t0)
        (3) inverse scaling, eta = eta0 / pow(t, power_t)
    eta0 : double
        The initial learning rate.
    power_t : double
        The exponent for inverse scaling learning rate.
    t : double
        Initial state of the learning rate. This value is equal to the
        iteration count except when the learning rate is set to `optimal`.
        Default: 1.0.
    intercept_decay : double
        An additional scaling factor for the intercept update.

    Returns
    -------
    weights : array, shape=[n_features]
        The fitted weight vector.
    intercept : float
        The fitted intercept term.

    """
    cdef WeightVector w = WeightVector(weights, intercept, fit_intercept,
                                       intercept_decay)

    # get the data information into easy vars
    cdef Py_ssize_t n_samples = dataset.n_samples
    cdef Py_ssize_t n_features = w.n_features

    # ``K==1`` for both regression and binary classification;
    # ``K==n_classes`` for multi-class classification.
    cdef int K = w.K

    # function pointer to concrete weight update
    cdef double (*weight_update)(LossFunction, WeightVector, DOUBLE *,
                                 INTEGER *, int, double, double,
                                 double *, double *, double, int, int)
    if K == 1:
        weight_update = &weight_update_binary
    else:
        weight_update = &weight_update_multinomial

    cdef DOUBLE *x_data_ptr = NULL
    cdef INTEGER *x_ind_ptr = NULL

    cdef double *class_weight_data_ptr = <double *>class_weight.data

    # helper variable
    cdef int xnnz
    cdef double eta = eta0
    cdef double sumloss = 0.0
    cdef double y = 0.0
    cdef double sample_weight
    cdef unsigned int count = 0
    cdef unsigned int epoch = 0
    cdef unsigned int i = 0

    # predictions array
    cdef double *p = NULL
    p = <double *> malloc(K * sizeof(double))
    for i in range(K):
        p[i] = 0.0

    # q array is only used for L1 regularization
    # it stores the cumulative L1 penalty for each feature
    cdef double *q = NULL
    if penalty_type == L1 or penalty_type == ELASTICNET:
        q = <double *> malloc(n_features * sizeof(double))
        for i in range(n_features):
            q[i] = 0.0
    # total cumulative L1 penalty
    cdef double u = 0.0

    if penalty_type == L2:
        rho = 1.0
    elif penalty_type == L1:
        rho = 0.0

    t_start = time()
    for epoch in range(n_iter):
        if verbose > 0:
            print("-- Epoch %d" % (epoch + 1))
        if shuffle:
            dataset.shuffle(seed)
        for i in range(n_samples):
            dataset.next(&x_data_ptr, &x_ind_ptr, &xnnz, &y,
                         &sample_weight)

            if learning_rate == OPTIMAL:
                eta = 1.0 / (alpha * t)
            elif learning_rate == INVSCALING:
                eta = eta0 / pow(t, power_t)

            weight_update(loss, w, x_data_ptr, x_ind_ptr, xnnz, y,
                          sample_weight, class_weight_data_ptr,
                          p, eta, fit_intercept, K)

            if penalty_type >= L2:
                w.scale(1.0 - (rho * eta * alpha))

            if penalty_type == L1 or penalty_type == ELASTICNET:
                u += ((1.0 - rho) * eta * alpha)
                l1penalty(w, q, x_ind_ptr, xnnz, u)
            t += 1
            count += 1

        # report epoch information
        if verbose > 0:
            print("Norm: %.2f, NNZs: %d, "\
            "T: %d, Avg. loss: %.6f" % (w.norm(),
                                        weights.nonzero()[0].shape[0],
                                        count,
                                        sumloss / count))
            print("Total training time: %.2f seconds." % (time() - t_start))

        # floating-point under-/overflow check.
        ## if np.any(np.isinf(weights)) or np.any(np.isnan(weights)) \
        ##    or np.any(np.isnan(intercept)) or np.any(np.isinf(intercept)):
        ##     raise ValueError("floating-point under-/overflow occured.")

    w.reset_wscale()

    free(p)
    if q != NULL:
        free(q)

    return weights, intercept


cdef double weight_update_binary(LossFunction loss, WeightVector w,
                                 DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                                 int xnnz, double y, double sample_weight,
                                 double *class_weights, double *p,
                                 double eta, int fit_intercept, int K):

    cdef double update = 0.0
    cdef double class_weight = class_weights[0]
    p[0] = w.dot(x_data_ptr, x_ind_ptr, xnnz, 0)

    if y > 0.0:
        class_weight = class_weights[1]

    update = eta * loss.dloss(p[0], y) * class_weight * sample_weight
    if update != 0.0:
        w.add(x_data_ptr, x_ind_ptr, xnnz, 0, -update)

    return 0.0  ##loss.loss(p[0], y)


cdef double weight_update_multinomial(LossFunction loss, WeightVector w,
                                      DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                                      int xnnz, double y, double sample_weight,
                                      double *class_weights,
                                      double *p, double eta, int fit_intercept,
                                      int K):
    cdef double update = 0.0
    cdef int k
    cdef double p_sum = 0.0

    for k in range(K):
        p[k] = exp(w.dot(x_data_ptr, x_ind_ptr, xnnz, k))
        p_sum += p[k]

    w.add(x_data_ptr, x_ind_ptr, xnnz, <int>y, eta)
    for k in range(K):
        update = -1.0 * eta * (p[k] / p_sum)
        w.add(x_data_ptr, x_ind_ptr, xnnz, k, update)

    return 0.0  ##log(p[<int>y] / p_sum)


cdef inline double max(double a, double b):
    return a if a >= b else b


cdef inline double min(double a, double b):
    return a if a <= b else b


cdef void l1penalty(WeightVector w, double *q,
                    INTEGER *x_ind_ptr, int xnnz, double u):
    """Apply the L1 penalty to each updated feature

    This implements the truncated gradient approach by
    [Tsuruoka, Y., Tsujii, J., and Ananiadou, S., 2009].
    """
    cdef double z = 0.0
    cdef int j = 0
    cdef int idx = 0
    cdef int n_features = w.n_features
    cdef int k = 0
    cdef int K = w.K
    cdef double wscale = w.wscale
    cdef double* w_data_ptr = w.w_data_ptr

    for j in range(xnnz):
        for k in range(K):
            idx = x_ind_ptr[j] + (k * n_features)
            z = w_data_ptr[idx]
            if (wscale * w_data_ptr[idx]) > 0.0:
                w_data_ptr[idx] = max(
                    0.0, w_data_ptr[idx] - ((u + q[idx]) / wscale))

            elif (wscale * w_data_ptr[idx]) < 0.0:
                w_data_ptr[idx] = min(
                    0.0, w_data_ptr[idx] + ((u - q[idx]) / wscale))

            q[idx] += (wscale * (w_data_ptr[idx] - z))
