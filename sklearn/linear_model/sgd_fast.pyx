# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Mathieu Blondel (partial_fit support)
#         Rob Zinkov (passive-aggressive)
#         Lars Buitinck
#
# License: BSD 3 clause


import numpy as np
import sys
from time import time

cimport cython
from libc.math cimport exp, log, sqrt, pow, fabs
cimport numpy as np
from numpy.math cimport INFINITY
cdef extern from "sgd_fast_helpers.h":
    bint skl_isfinite(double) nogil

from sklearn.utils.weight_vector cimport WeightVector
from sklearn.utils.seq_dataset cimport SequentialDataset

np.import_array()

# Penalty constants
DEF NO_PENALTY = 0
DEF L1 = 1
DEF L2 = 2
DEF ELASTICNET = 3

# Learning rate constants
DEF CONSTANT = 1
DEF OPTIMAL = 2
DEF INVSCALING = 3
DEF ADAPTIVE = 4
DEF PA1 = 5
DEF PA2 = 6



# ----------------------------------------
# Extension Types for Loss Functions
# ----------------------------------------

cdef class LossFunction:
    """Base class for convex loss functions"""

    cdef double loss(self, double p, double y) nogil:
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
        return 0.

    def dloss(self, double p, double y):
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
            The derivative of the loss function with regards to `p`.
        """
        return self._dloss(p, y)

    cdef double _dloss(self, double p, double y) nogil:
        # Implementation of dloss; separate function because cpdef and nogil
        # can't be combined.
        return 0.


cdef class Regression(LossFunction):
    """Base class for loss functions for regression"""

    cdef double loss(self, double p, double y) nogil:
        return 0.

    cdef double _dloss(self, double p, double y) nogil:
        return 0.


cdef class Classification(LossFunction):
    """Base class for loss functions for classification"""

    cdef double loss(self, double p, double y) nogil:
        return 0.

    cdef double _dloss(self, double p, double y) nogil:
        return 0.


cdef class ModifiedHuber(Classification):
    """Modified Huber loss for binary classification with y in {-1, 1}

    This is equivalent to quadratically smoothed SVM with gamma = 2.

    See T. Zhang 'Solving Large Scale Linear Prediction Problems Using
    Stochastic Gradient Descent', ICML'04.
    """
    cdef double loss(self, double p, double y) nogil:
        cdef double z = p * y
        if z >= 1.0:
            return 0.0
        elif z >= -1.0:
            return (1.0 - z) * (1.0 - z)
        else:
            return -4.0 * z

    cdef double _dloss(self, double p, double y) nogil:
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

    cdef double loss(self, double p, double y) nogil:
        cdef double z = p * y
        if z <= self.threshold:
            return self.threshold - z
        return 0.0

    cdef double _dloss(self, double p, double y) nogil:
        cdef double z = p * y
        if z <= self.threshold:
            return -y
        return 0.0

    def __reduce__(self):
        return Hinge, (self.threshold,)


cdef class SquaredHinge(Classification):
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

    cdef double loss(self, double p, double y) nogil:
        cdef double z = self.threshold - p * y
        if z > 0:
            return z * z
        return 0.0

    cdef double _dloss(self, double p, double y) nogil:
        cdef double z = self.threshold - p * y
        if z > 0:
            return -2 * y * z
        return 0.0

    def __reduce__(self):
        return SquaredHinge, (self.threshold,)


cdef class Log(Classification):
    """Logistic regression loss for binary classification with y in {-1, 1}"""

    cdef double loss(self, double p, double y) nogil:
        cdef double z = p * y
        # approximately equal and saves the computation of the log
        if z > 18:
            return exp(-z)
        if z < -18:
            return -z
        return log(1.0 + exp(-z))

    cdef double _dloss(self, double p, double y) nogil:
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
    cdef double loss(self, double p, double y) nogil:
        return 0.5 * (p - y) * (p - y)

    cdef double _dloss(self, double p, double y) nogil:
        return p - y

    def __reduce__(self):
        return SquaredLoss, ()


cdef class Huber(Regression):
    """Huber regression loss

    Variant of the SquaredLoss that is robust to outliers (quadratic near zero,
    linear in for large errors).

    https://en.wikipedia.org/wiki/Huber_Loss_Function
    """

    cdef double c

    def __init__(self, double c):
        self.c = c

    cdef double loss(self, double p, double y) nogil:
        cdef double r = p - y
        cdef double abs_r = fabs(r)
        if abs_r <= self.c:
            return 0.5 * r * r
        else:
            return self.c * abs_r - (0.5 * self.c * self.c)

    cdef double _dloss(self, double p, double y) nogil:
        cdef double r = p - y
        cdef double abs_r = fabs(r)
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

    cdef double loss(self, double p, double y) nogil:
        cdef double ret = fabs(y - p) - self.epsilon
        return ret if ret > 0 else 0

    cdef double _dloss(self, double p, double y) nogil:
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

    cdef double loss(self, double p, double y) nogil:
        cdef double ret = fabs(y - p) - self.epsilon
        return ret * ret if ret > 0 else 0

    cdef double _dloss(self, double p, double y) nogil:
        cdef double z
        z = y - p
        if z > self.epsilon:
            return -2 * (z - self.epsilon)
        elif z < -self.epsilon:
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
              double l1_ratio,
              SequentialDataset dataset,
              np.ndarray[unsigned char, ndim=1, mode='c'] validation_mask,
              bint early_stopping, estimator,
              int n_iter_no_change,
              int max_iter, double tol, int fit_intercept,
              int verbose, bint shuffle, np.uint32_t seed,
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
    C : float
        Maximum step size for passive aggressive.
    l1_ratio : float
        The Elastic Net mixing parameter, with 0 <= l1_ratio <= 1.
        l1_ratio=0 corresponds to L2 penalty, l1_ratio=1 to L1.
    dataset : SequentialDataset
        A concrete ``SequentialDataset`` object.
    validation_mask : ndarray[unsigned char, ndim=1]
        Equal to True on the validation set.
    early_stopping : boolean
        Whether to use a stopping criterion based on the validation set.
    estimator : BaseSGD
        A concrete object inheriting from ``BaseSGD``.
        Used only if early_stopping is True.
    n_iter_no_change : int
        Number of iteration with no improvement to wait before stopping.
    max_iter : int
        The maximum number of iterations (epochs).
    tol: double
        The tolerance for the stopping criterion.
    fit_intercept : int
        Whether or not to fit the intercept (1 or 0).
    verbose : int
        Print verbose output; 0 for quite.
    shuffle : boolean
        Whether to shuffle the training data before each epoch.
    weight_pos : float
        The weight of the positive class.
    weight_neg : float
        The weight of the negative class.
    seed : np.uint32_t
        Seed of the pseudorandom number generator used to shuffle the data.
    learning_rate : int
        The learning rate:
        (1) constant, eta = eta0
        (2) optimal, eta = 1.0/(alpha * t).
        (3) inverse scaling, eta = eta0 / pow(t, power_t)
        (4) adaptive decrease
        (5) Passive Aggressive-I, eta = min(alpha, loss/norm(x))
        (6) Passive Aggressive-II, eta = 1.0 / (norm(x) + 0.5*alpha)
    eta0 : double
        The initial learning rate.
    power_t : double
        The exponent for inverse scaling learning rate.
    t : double
        Initial state of the learning rate. This value is equal to the
        iteration count except when the learning rate is set to `optimal`.
        Default: 1.0.
    intercept_decay : double
        The decay ratio of intercept, used in updating intercept.

    Returns
    -------
    weights : array, shape=[n_features]
        The fitted weight vector.
    intercept : float
        The fitted intercept term.
    n_iter_ : int
        The actual number of iter (epochs).
    """
    standard_weights, standard_intercept,\
        _, _, n_iter_ = _plain_sgd(weights,
                                   intercept,
                                   None,
                                   0,
                                   loss,
                                   penalty_type,
                                   alpha, C,
                                   l1_ratio,
                                   dataset,
                                   validation_mask,
                                   early_stopping,
                                   estimator,
                                   n_iter_no_change,
                                   max_iter, tol, fit_intercept,
                                   verbose, shuffle, seed,
                                   weight_pos, weight_neg,
                                   learning_rate, eta0,
                                   power_t,
                                   t,
                                   intercept_decay,
                                   0)
    return standard_weights, standard_intercept, n_iter_


def average_sgd(np.ndarray[double, ndim=1, mode='c'] weights,
                double intercept,
                np.ndarray[double, ndim=1, mode='c'] average_weights,
                double average_intercept,
                LossFunction loss,
                int penalty_type,
                double alpha, double C,
                double l1_ratio,
                SequentialDataset dataset,
                np.ndarray[unsigned char, ndim=1, mode='c'] validation_mask,
                bint early_stopping, estimator,
                int n_iter_no_change,
                int max_iter, double tol, int fit_intercept,
                int verbose, bint shuffle, np.uint32_t seed,
                double weight_pos, double weight_neg,
                int learning_rate, double eta0,
                double power_t,
                double t=1.0,
                double intercept_decay=1.0,
                int average=1):
    """Average SGD for generic loss functions and penalties.

    Parameters
    ----------
    weights : ndarray[double, ndim=1]
        The allocated coef_ vector.
    intercept : double
        The initial intercept.
    average_weights : ndarray[double, ndim=1]
        The average weights as computed for ASGD
    average_intercept : double
        The average intercept for ASGD
    loss : LossFunction
        A concrete ``LossFunction`` object.
    penalty_type : int
        The penalty 2 for L2, 1 for L1, and 3 for Elastic-Net.
    alpha : float
        The regularization parameter.
    C : float
        Maximum step size for passive aggressive.
    l1_ratio : float
        The Elastic Net mixing parameter, with 0 <= l1_ratio <= 1.
        l1_ratio=0 corresponds to L2 penalty, l1_ratio=1 to L1.
    dataset : SequentialDataset
        A concrete ``SequentialDataset`` object.
    validation_mask : ndarray[unsigned char, ndim=1]
        Equal to True on the validation set.
    early_stopping : boolean
        Whether to use a stopping criterion based on the validation set.
    estimator : BaseSGD
        A concrete object inheriting from ``BaseSGD``.
        Used only if early_stopping is True.
    n_iter_no_change : int
        Number of iteration with no improvement to wait before stopping.
    max_iter : int
        The maximum number of iterations (epochs).
    tol: double
        The tolerance for the stopping criterion.
    fit_intercept : int
        Whether or not to fit the intercept (1 or 0).
    verbose : int
        Print verbose output; 0 for quite.
    shuffle : boolean
        Whether to shuffle the training data before each epoch.
    weight_pos : float
        The weight of the positive class.
    weight_neg : float
        The weight of the negative class.
    seed : np.uint32_t
        Seed of the pseudorandom number generator used to shuffle the data.
    learning_rate : int
        The learning rate:
        (1) constant, eta = eta0
        (2) optimal, eta = 1.0/(alpha * t).
        (3) inverse scaling, eta = eta0 / pow(t, power_t)
        (4) adaptive decrease
        (5) Passive Aggressive-I, eta = min(alpha, loss/norm(x))
        (6) Passive Aggressive-II, eta = 1.0 / (norm(x) + 0.5*alpha)
    eta0 : double
        The initial learning rate.
    power_t : double
        The exponent for inverse scaling learning rate.
    t : double
        Initial state of the learning rate. This value is equal to the
        iteration count except when the learning rate is set to `optimal`.
        Default: 1.0.
    average : int
        The number of iterations before averaging starts. average=1 is
        equivalent to averaging for all iterations.

    Returns
    -------
    weights : array, shape=[n_features]
        The fitted weight vector.
    intercept : float
        The fitted intercept term.
    average_weights : array shape=[n_features]
        The averaged weights across iterations
    average_intercept : float
        The averaged intercept across iterations
    n_iter_ : int
        The actual number of iter (epochs).
    """
    return _plain_sgd(weights,
                      intercept,
                      average_weights,
                      average_intercept,
                      loss,
                      penalty_type,
                      alpha, C,
                      l1_ratio,
                      dataset,
                      validation_mask,
                      early_stopping,
                      estimator,
                      n_iter_no_change,
                      max_iter, tol, fit_intercept,
                      verbose, shuffle, seed,
                      weight_pos, weight_neg,
                      learning_rate, eta0,
                      power_t,
                      t,
                      intercept_decay,
                      average)


def _plain_sgd(np.ndarray[double, ndim=1, mode='c'] weights,
               double intercept,
               np.ndarray[double, ndim=1, mode='c'] average_weights,
               double average_intercept,
               LossFunction loss,
               int penalty_type,
               double alpha, double C,
               double l1_ratio,
               SequentialDataset dataset,
               np.ndarray[unsigned char, ndim=1, mode='c'] validation_mask,
               bint early_stopping, estimator,
               int n_iter_no_change,
               int max_iter, double tol, int fit_intercept,
               int verbose, bint shuffle, np.uint32_t seed,
               double weight_pos, double weight_neg,
               int learning_rate, double eta0,
               double power_t,
               double t=1.0,
               double intercept_decay=1.0,
               int average=0):

    # get the data information into easy vars
    cdef Py_ssize_t n_samples = dataset.n_samples
    cdef Py_ssize_t n_features = weights.shape[0]

    cdef WeightVector w = WeightVector(weights, average_weights)
    cdef double* w_ptr = &weights[0]
    cdef double *x_data_ptr = NULL
    cdef int *x_ind_ptr = NULL
    cdef double* ps_ptr = NULL

    # helper variables
    cdef int no_improvement_count = 0
    cdef bint infinity = False
    cdef int xnnz
    cdef double eta = 0.0
    cdef double p = 0.0
    cdef double update = 0.0
    cdef double sumloss = 0.0
    cdef double score = 0.0
    cdef double best_loss = INFINITY
    cdef double best_score = -INFINITY
    cdef double y = 0.0
    cdef double sample_weight
    cdef double class_weight = 1.0
    cdef unsigned int count = 0
    cdef unsigned int epoch = 0
    cdef unsigned int i = 0
    cdef int is_hinge = isinstance(loss, Hinge)
    cdef double optimal_init = 0.0
    cdef double dloss = 0.0
    cdef double MAX_DLOSS = 1e12
    cdef double max_change = 0.0
    cdef double max_weight = 0.0

    cdef long long sample_index
    cdef unsigned char [:] validation_mask_view = validation_mask

    # q vector is only used for L1 regularization
    cdef np.ndarray[double, ndim = 1, mode = "c"] q = None
    cdef double * q_data_ptr = NULL
    if penalty_type == L1 or penalty_type == ELASTICNET:
        q = np.zeros((n_features,), dtype=np.float64, order="c")
        q_data_ptr = <double * > q.data
    cdef double u = 0.0

    if penalty_type == L2:
        l1_ratio = 0.0
    elif penalty_type == L1:
        l1_ratio = 1.0

    eta = eta0

    if learning_rate == OPTIMAL:
        typw = np.sqrt(1.0 / np.sqrt(alpha))
        # computing eta0, the initial learning rate
        initial_eta0 = typw / max(1.0, loss.dloss(-typw, 1.0))
        # initialize t such that eta at first sample equals eta0
        optimal_init = 1.0 / (initial_eta0 * alpha)

    t_start = time()
    with nogil:
        for epoch in range(max_iter):
            sumloss = 0
            if verbose > 0:
                with gil:
                    print("-- Epoch %d" % (epoch + 1))
            if shuffle:
                dataset.shuffle(seed)
            for i in range(n_samples):
                dataset.next(&x_data_ptr, &x_ind_ptr, &xnnz,
                             &y, &sample_weight)

                sample_index = dataset.index_data_ptr[dataset.current_index]
                if validation_mask_view[sample_index]:
                    # do not learn on the validation set
                    continue

                p = w.dot(x_data_ptr, x_ind_ptr, xnnz) + intercept
                if learning_rate == OPTIMAL:
                    eta = 1.0 / (alpha * (optimal_init + t - 1))
                elif learning_rate == INVSCALING:
                    eta = eta0 / pow(t, power_t)

                if verbose or not early_stopping:
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
                    dloss = loss._dloss(p, y)
                    # clip dloss with large values to avoid numerical
                    # instabilities
                    if dloss < -MAX_DLOSS:
                        dloss = -MAX_DLOSS
                    elif dloss > MAX_DLOSS:
                        dloss = MAX_DLOSS
                    update = -eta * dloss

                if learning_rate >= PA1:
                    if is_hinge:
                        # classification
                        update *= y
                    elif y - p < 0:
                        # regression
                        update *= -1

                update *= class_weight * sample_weight

                if penalty_type >= L2:
                    # do not scale to negative values when eta or alpha are too
                    # big: instead set the weights to zero
                    w.scale(max(0, 1.0 - ((1.0 - l1_ratio) * eta * alpha)))
                if update != 0.0:
                    w.add(x_data_ptr, x_ind_ptr, xnnz, update)
                    if fit_intercept == 1:
                        intercept += update * intercept_decay

                if 0 < average <= t:
                    # compute the average for the intercept and update the
                    # average weights, this is done regardless as to whether
                    # the update is 0

                    w.add_average(x_data_ptr, x_ind_ptr, xnnz,
                                  update, (t - average + 1))
                    average_intercept += ((intercept - average_intercept) /
                                          (t - average + 1))

                if penalty_type == L1 or penalty_type == ELASTICNET:
                    u += (l1_ratio * eta * alpha)
                    l1penalty(w, q_data_ptr, x_ind_ptr, xnnz, u)

                t += 1
                count += 1

            # report epoch information
            if verbose > 0:
                with gil:
                    print("Norm: %.2f, NNZs: %d, Bias: %.6f, T: %d, "
                          "Avg. loss: %f"
                          % (w.norm(), weights.nonzero()[0].shape[0],
                             intercept, count, sumloss / n_samples))
                    print("Total training time: %.2f seconds."
                          % (time() - t_start))

            # floating-point under-/overflow check.
            if (not skl_isfinite(intercept)
                or any_nonfinite(<double *>weights.data, n_features)):
                infinity = True
                break

            # evaluate the score on the validation set
            if early_stopping:
                with gil:
                    score = estimator._validation_score(weights, intercept)
                if tol > -INFINITY and score < best_score + tol:
                    no_improvement_count += 1
                else:
                    no_improvement_count = 0
                if score > best_score:
                    best_score = score
            # or evaluate the loss on the training set
            else:
                if tol > -INFINITY and sumloss > best_loss - tol * n_samples:
                    no_improvement_count += 1
                else:
                    no_improvement_count = 0
                if sumloss < best_loss:
                    best_loss = sumloss

            # if there is no improvement several times in a row
            if no_improvement_count >= n_iter_no_change:
                if learning_rate == ADAPTIVE and eta > 1e-6:
                    eta = eta / 5
                    no_improvement_count = 0
                else:
                    if verbose:
                        with gil:
                            print("Convergence after %d epochs took %.2f "
                                  "seconds" % (epoch + 1, time() - t_start))
                    break

    if infinity:
        raise ValueError(("Floating-point under-/overflow occurred at epoch"
                          " #%d. Scaling input data with StandardScaler or"
                          " MinMaxScaler might help.") % (epoch + 1))

    w.reset_wscale()

    return weights, intercept, average_weights, average_intercept, epoch + 1


cdef bint any_nonfinite(double *w, int n) nogil:
    for i in range(n):
        if not skl_isfinite(w[i]):
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
