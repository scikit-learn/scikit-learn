# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# Licence: BSD 3 clause


import numpy as np
import sys
from time import time

cimport cython
from libc.math cimport exp, log, sqrt, pow, fabs
cimport numpy as np
cdef extern from "numpy/npy_math.h":
    bint isfinite "npy_isfinite"(double) nogil

from sklearn.utils.weight_vector cimport WeightVector
from sklearn.utils.seq_dataset cimport SequentialDataset

np.import_array()


# Penalty constans
DEF NO_PENALTY = 0
DEF L1 = 1
DEF L2 = 2
DEF ELASTICNET = 3

# Learning rate constants
DEF CONSTANT = 1
DEF OPTIMAL = 2
DEF INVSCALING = 3
DEF PA1 = 4
DEF PA2 = 5

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


cdef class SquaredHinge(LossFunction):
    """Squared Hinge loss for binary classification tasks with y in {-1,1}

    Parameters
    ----------

    threshold : float > 0.0
        Margin threshold. When threshold=1.0, one gets the loss used by
        (quadratically penalized) SVM.
    """

    cdef double threshold

    def __init__(self, double threshold=1.0):
        self.threshold = threshold

    cpdef double loss(self, double p, double y):
        cdef double z = self.threshold - p * y
        if z > 0:
            return z * z
        return 0.0

    cpdef double dloss(self, double p, double y):
        cdef double z = self.threshold - p * y
        if z > 0:
            return -2 * y * z
        return 0.0

    def __reduce__(self):
        return SquaredHinge, (self.threshold,)


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


cdef class EpsilonInsensitive(Regression):
    """Epsilon-Insensitive loss (used by SVR).

    loss = max(0, |y - p| - epsilon)
    """

    cdef double epsilon

    def __init__(self, double epsilon):
        self.epsilon = epsilon

    cpdef double loss(self, double p, double y):
        cdef double ret = abs(y - p) - self.epsilon
        return ret if ret > 0 else 0

    cpdef double dloss(self, double p, double y):
        if y - p > self.epsilon:
            return -1
        elif p - y > self.epsilon:
            return 1
        else:
            return 0

    def __reduce__(self):
        return EpsilonInsensitive, (self.epsilon,)


cdef class SquaredEpsilonInsensitive(Regression):
    """Epsilon-Insensitive loss.

    loss = max(0, |y - p| - epsilon)^2
    """

    cdef double epsilon

    def __init__(self, double epsilon):
        self.epsilon = epsilon

    cpdef double loss(self, double p, double y):
        cdef double ret = abs(y - p) - self.epsilon
        return ret * ret if ret > 0 else 0

    cpdef double dloss(self, double p, double y):
        cdef double z
        z = y - p
        if z > self.epsilon:
            return -2 * (z - self.epsilon)
        elif z < self.epsilon:
            return 2 * (-z - self.epsilon)
        else:
            return 0

    def __reduce__(self):
        return SquaredEpsilonInsensitive, (self.epsilon,)


def plain_sgd(np.ndarray[double, ndim=1, mode='c'] weights,
              double intercept,
              LossFunction loss,
              int penalty_type,
              double alpha, double C,
              double rho,
              SequentialDataset dataset,
              int n_iter, int fit_intercept,
              int verbose, int shuffle, seed,
              double weight_pos, double weight_neg,
              int learning_rate, double eta0,
              double power_t,
              double t=1.0,
              double intercept_decay=1.0):
    """Plain SGD for generic loss functions and penalties.

    Parameters
    ----------
    weights : ndarray[double, ndim=1]
        The allocated coef_ vector.
    intercept : double
        The initial intercept.
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
    weight_pos : float
        The weight of the positive class.
    weight_neg : float
        The weight of the negative class.
    seed : int or RandomState object
        The seed of the pseudo random number generator to use when
        shuffling the data
    learning_rate : int
        The learning rate:
        (1) constant, eta = eta0
        (2) optimal, eta = 1.0/(t+t0)
        (3) inverse scaling, eta = eta0 / pow(t, power_t)
        (4) Passive Agressive-I, eta = min(alpha, loss/norm(x))
        (5) Passive Agressive-II, eta = 1.0 / (norm(x) + 0.5*alpha)
    eta0 : double
        The initial learning rate.
    power_t : double
        The exponent for inverse scaling learning rate.
    t : double
        Initial state of the learning rate. This value is equal to the
        iteration count except when the learning rate is set to `optimal`.
        Default: 1.0.

    Returns
    -------
    weights : array, shape=[n_features]
        The fitted weight vector.
    intercept : float
        The fitted intercept term.

    """

    # get the data information into easy vars
    cdef Py_ssize_t n_samples = dataset.n_samples
    cdef Py_ssize_t n_features = weights.shape[0]

    cdef WeightVector w = WeightVector(weights)

    cdef double *x_data_ptr = NULL
    cdef int *x_ind_ptr = NULL

    # helper variable
    cdef int xnnz
    cdef double eta = 0.0
    cdef double p = 0.0
    cdef double update = 0.0
    cdef double sumloss = 0.0
    cdef double y = 0.0
    cdef double sample_weight
    cdef double class_weight = 1.0
    cdef unsigned int count = 0
    cdef unsigned int epoch = 0
    cdef unsigned int i = 0
    cdef int is_hinge = isinstance(loss, Hinge)

    # q vector is only used for L1 regularization
    cdef np.ndarray[double, ndim = 1, mode = "c"] q = None
    cdef double * q_data_ptr = NULL
    if penalty_type == L1 or penalty_type == ELASTICNET:
        q = np.zeros((n_features,), dtype=np.float64, order="c")
        q_data_ptr = <double * > q.data
    cdef double u = 0.0

    if penalty_type == L2:
        rho = 1.0
    elif penalty_type == L1:
        rho = 0.0

    eta = eta0

    t_start = time()
    for epoch in range(n_iter):
        if verbose > 0:
            print("-- Epoch %d" % (epoch + 1))
        if shuffle:
            dataset.shuffle(seed)
        for i in range(n_samples):
            dataset.next( & x_data_ptr, & x_ind_ptr, & xnnz, & y,
                         & sample_weight)

            p = w.dot(x_data_ptr, x_ind_ptr, xnnz) + intercept

            if learning_rate == OPTIMAL:
                eta = 1.0 / (alpha * t)
            elif learning_rate == INVSCALING:
                eta = eta0 / pow(t, power_t)

            if verbose > 0:
                sumloss += loss.loss(p, y)

            if y > 0.0:
                class_weight = weight_pos
            else:
                class_weight = weight_neg

            if learning_rate == PA1:
                update = sqnorm(x_data_ptr, x_ind_ptr, xnnz)
                if update == 0:
                    continue
                update = min(C, loss.loss(p, y) / update)
            elif learning_rate == PA2:
                update = sqnorm(x_data_ptr, x_ind_ptr, xnnz)
                update = loss.loss(p, y) / (update + 0.5 / C)
            else:
                update = -eta * loss.dloss(p, y)

            if learning_rate >= PA1:
                if is_hinge:
                    # classification
                    update *= y
                elif y - p < 0:
                    # regression
                    update *= -1

            update *= class_weight * sample_weight

            if penalty_type >= L2:
                w.scale(1.0 - (rho * eta * alpha))
            if update != 0.0:
                w.add(x_data_ptr, x_ind_ptr, xnnz, update)
                if fit_intercept == 1:
                    intercept += update * intercept_decay

            if penalty_type == L1 or penalty_type == ELASTICNET:
                u += ((1.0 - rho) * eta * alpha)
                l1penalty(w, q_data_ptr, x_ind_ptr, xnnz, u)
            t += 1
            count += 1

        # report epoch information
        if verbose > 0:
            print("Norm: %.2f, NNZs: %d, "
            "Bias: %.6f, T: %d, Avg. loss: %.6f" % (w.norm(),
                                                    weights.nonzero()[0].shape[0],
                                                    intercept, count,
                                                    sumloss / count))
            print("Total training time: %.2f seconds." % (time() - t_start))

        # floating-point under-/overflow check.
        if (not isfinite(intercept)
            or any_nonfinite(<double *>weights.data, n_features)):
            raise ValueError("floating-point under-/overflow occurred.")

    w.reset_wscale()

    return weights, intercept


cdef bint any_nonfinite(double *w, int n):
    cdef int i

    for i in range(n):
        if not isfinite(w[i]):
            return True
    return 0


cdef double sqnorm(double * x_data_ptr, int * x_ind_ptr, int xnnz) nogil:
    cdef double x_norm = 0.0
    cdef int j
    cdef double z
    for j in range(xnnz):
        z = x_data_ptr[j]
        x_norm += z * z
    return x_norm


cdef void l1penalty(WeightVector w, double * q_data_ptr,
                    int *x_ind_ptr, int xnnz, double u) nogil:
    """Apply the L1 penalty to each updated feature

    This implements the truncated gradient approach by
    [Tsuruoka, Y., Tsujii, J., and Ananiadou, S., 2009].
    """
    cdef double z = 0.0
    cdef int j = 0
    cdef int idx = 0
    cdef double wscale = w.wscale
    cdef double *w_data_ptr = w.w_data_ptr
    for j in range(xnnz):
        idx = x_ind_ptr[j]
        z = w_data_ptr[idx]
        if wscale * w_data_ptr[idx] > 0.0:
            w_data_ptr[idx] = max(
                0.0, w_data_ptr[idx] - ((u + q_data_ptr[idx]) / wscale))

        elif wscale * w_data_ptr[idx] < 0.0:
            w_data_ptr[idx] = min(
                0.0, w_data_ptr[idx] + ((u - q_data_ptr[idx]) / wscale))

        q_data_ptr[idx] += wscale * (w_data_ptr[idx] - z)
