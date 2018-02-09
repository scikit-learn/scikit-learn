# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Authors: Danny Sullivan <dbsullivan23@gmail.com>
#          Tom Dupre la Tour <tom.dupre-la-tour@m4x.org>
#          Arthur Mensch <arthur.mensch@m4x.org
#
# License: BSD 3 clause
cimport numpy as np
import numpy as np
from libc.math cimport fabs, exp, log
from libc.time cimport time, time_t

from .sgd_fast cimport LossFunction
from .sgd_fast cimport Log, SquaredLoss
from ..utils.seq_dataset cimport SequentialDataset

from libc.stdio cimport printf

cdef extern from "sgd_fast_helpers.h":
    bint skl_isfinite(double) nogil


cdef inline double fmax(double x, double y) nogil:
    if x > y:
        return x
    return y


cdef double _logsumexp(double* arr, int n_classes) nogil:
    """Computes the sum of arr assuming arr is in the log domain.

    Returns log(sum(exp(arr))) while minimizing the possibility of
    over/underflow.
    """
    # Use the max to normalize, as with the log this is what accumulates
    # the less errors
    cdef double vmax = arr[0]
    cdef double out = 0.0
    cdef int i

    for i in range(1, n_classes):
        if vmax < arr[i]:
            vmax = arr[i]

    for i in range(n_classes):
        out += exp(arr[i] - vmax)

    return log(out) + vmax


cdef class MultinomialLogLoss:
    cdef double _loss(self, double* prediction, double y, int n_classes,
                      double sample_weight) nogil:
        r"""Multinomial Logistic regression loss.

        The multinomial logistic loss for one sample is:
        loss = - sw \sum_c \delta_{y,c} (prediction[c] - logsumexp(prediction))
             = sw (logsumexp(prediction) - prediction[y])

        where:
            prediction = dot(x_sample, weights) + intercept
            \delta_{y,c} = 1 if (y == c) else 0
            sw = sample_weight

        Parameters
        ----------
        prediction : pointer to a np.ndarray[double] of shape (n_classes,)
            Prediction of the multinomial classifier, for current sample.

        y : double, between 0 and n_classes - 1
            Indice of the correct class for current sample (i.e. label encoded).

        n_classes : integer
            Total number of classes.

        sample_weight : double
            Weight of current sample.

        Returns
        -------
        loss : double
            Multinomial loss for current sample.

        Reference
        ---------
        Bishop, C. M. (2006). Pattern recognition and machine learning.
        Springer. (Chapter 4.3.4)
        """
        cdef double logsumexp_prediction = _logsumexp(prediction, n_classes)
        cdef double loss

        # y is the indice of the correct class of current sample.
        loss = (logsumexp_prediction - prediction[int(y)]) * sample_weight
        return loss

    cdef void _dloss(self, double* prediction, double y, int n_classes,
                     double sample_weight, double* gradient_ptr) nogil:
        r"""Multinomial Logistic regression gradient of the loss.

        The gradient of the multinomial logistic loss with respect to a class c,
        and for one sample is:
        grad_c = - sw * (p[c] - \delta_{y,c})

        where:
            p[c] = exp(logsumexp(prediction) - prediction[c])
            prediction = dot(sample, weights) + intercept
            \delta_{y,c} = 1 if (y == c) else 0
            sw = sample_weight

        Note that to obtain the true gradient, this value has to be multiplied
        by the sample vector x.

        Parameters
        ----------
        prediction : pointer to a np.ndarray[double] of shape (n_classes,)
            Prediction of the multinomial classifier, for current sample.

        y : double, between 0 and n_classes - 1
            Indice of the correct class for current sample (i.e. label encoded)

        n_classes : integer
            Total number of classes.

        sample_weight : double
            Weight of current sample.

        gradient_ptr : pointer to a np.ndarray[double] of shape (n_classes,)
            Gradient vector to be filled.

        Reference
        ---------
        Bishop, C. M. (2006). Pattern recognition and machine learning.
        Springer. (Chapter 4.3.4)
        """
        cdef double logsumexp_prediction = _logsumexp(prediction, n_classes)
        cdef int class_ind

        for class_ind in range(n_classes):
            gradient_ptr[class_ind] = exp(prediction[class_ind] -
                                          logsumexp_prediction)

            # y is the indice of the correct class of current sample.
            if class_ind == y:
                gradient_ptr[class_ind] -= 1.0

            gradient_ptr[class_ind] *= sample_weight

    def __reduce__(self):
        return MultinomialLogLoss, ()


def _multinomial_grad_loss_all_samples(
        SequentialDataset dataset,
        np.ndarray[double, ndim=2, mode='c'] weights_array,
        np.ndarray[double, ndim=1, mode='c'] intercept_array,
        int n_samples, int n_features, int n_classes):
    """Compute multinomial gradient and loss across all samples.

    Used for testing purpose only.
    """
    cdef double* weights = <double * >weights_array.data
    cdef double* intercept = <double * >intercept_array.data

    cdef double *x_data_ptr = NULL
    cdef int *x_ind_ptr = NULL
    cdef int xnnz = -1
    cdef double y
    cdef double sample_weight

    cdef double wscale = 1.0
    cdef int i, j, class_ind, feature_ind
    cdef double val
    cdef double sum_loss = 0.0

    cdef MultinomialLogLoss multiloss = MultinomialLogLoss()

    cdef np.ndarray[double, ndim=2] sum_gradient_array = \
        np.zeros((n_features, n_classes), dtype=np.double, order="c")
    cdef double* sum_gradient = <double*> sum_gradient_array.data

    cdef np.ndarray[double, ndim=1] prediction_array = \
        np.zeros(n_classes, dtype=np.double, order="c")
    cdef double* prediction = <double*> prediction_array.data

    cdef np.ndarray[double, ndim=1] gradient_array = \
        np.zeros(n_classes, dtype=np.double, order="c")
    cdef double* gradient = <double*> gradient_array.data

    with nogil:
        for i in range(n_samples):
            # get next sample on the dataset
            dataset.next(&x_data_ptr, &x_ind_ptr, &xnnz,
                         &y, &sample_weight)

            # prediction of the multinomial classifier for the sample
            predict_sample(x_data_ptr, x_ind_ptr, xnnz, weights, wscale,
                           intercept, prediction, n_classes)

            # compute the gradient for this sample, given the prediction
            multiloss._dloss(prediction, y, n_classes, sample_weight, gradient)

            # compute the loss for this sample, given the prediction
            sum_loss += multiloss._loss(prediction, y, n_classes, sample_weight)

            # update the sum of the gradient
            for j in range(xnnz):
                feature_ind = x_ind_ptr[j]
                val = x_data_ptr[j]
                for class_ind in range(n_classes):
                    sum_gradient[feature_ind * n_classes + class_ind] += \
                        gradient[class_ind] * val

    return sum_loss, sum_gradient_array


cdef inline double _soft_thresholding(double x, double shrinkage) nogil:
    return fmax(x - shrinkage, 0) - fmax(- x - shrinkage, 0)


def sag(SequentialDataset dataset,
        np.ndarray[double, ndim=2, mode='c'] weights_array,
        np.ndarray[double, ndim=1, mode='c'] intercept_array,
        int n_samples,
        int n_features,
        int n_classes,
        double tol,
        int max_iter,
        str loss_function,
        double step_size,
        double alpha,
        double beta,
        np.ndarray[double, ndim=2, mode='c'] sum_gradient_init,
        np.ndarray[double, ndim=2, mode='c'] gradient_memory_init,
        np.ndarray[bint, ndim=1, mode='c'] seen_init,
        int num_seen,
        bint fit_intercept,
        np.ndarray[double, ndim=1, mode='c'] intercept_sum_gradient_init,
        double intercept_decay,
        bint saga,
        bint verbose):
    """Stochastic Average Gradient (SAG) and SAGA solvers.

    Used in Ridge and LogisticRegression.

    Reference
    ---------
    Schmidt, M., Roux, N. L., & Bach, F. (2013).
    Minimizing finite sums with the stochastic average gradient
    https://hal.inria.fr/hal-00860051/document
    (section 4.3)

    Defazio, A., Bach, F., Lacoste-Julien, S. (2014),
    SAGA: A Fast Incremental Gradient Method With Support
    for Non-Strongly Convex Composite Objectives
    https://arxiv.org/abs/1407.0202

    """
    # the data pointer for x, the current sample
    cdef double *x_data_ptr = NULL
    # the index pointer for the column of the data
    cdef int *x_ind_ptr = NULL
    # the number of non-zero features for current sample
    cdef int xnnz = -1
    # the label value for current sample
    cdef double y
    # the sample weight
    cdef double sample_weight

    # helper variable for indexes
    cdef int f_idx, s_idx, feature_ind, class_ind, j
    # the number of pass through all samples
    cdef int n_iter = 0
    # helper to track iterations through samples
    cdef int sample_itr
    # the index (row number) of the current sample
    cdef int sample_ind

    # the maximum change in weights, used to compute stopping criteria
    cdef double max_change
    # a holder variable for the max weight, used to compute stopping criteria
    cdef double max_weight

    # the start time of the fit
    cdef time_t start_time
    # the end time of the fit
    cdef time_t end_time

    # precomputation since the step size does not change in this implementation
    cdef double wscale_update = 1.0 - step_size * alpha

    # vector of booleans indicating whether this sample has been seen
    cdef bint* seen = <bint*> seen_init.data

    # helper for cumulative sum
    cdef double cum_sum

    # the pointer to the coef_ or weights
    cdef double* weights = <double * >weights_array.data
    # the pointer to the intercept_array
    cdef double* intercept = <double * >intercept_array.data

    # the pointer to the intercept_sum_gradient
    cdef double* intercept_sum_gradient = \
        <double * >intercept_sum_gradient_init.data

    # the sum of gradients for each feature
    cdef double* sum_gradient = <double*> sum_gradient_init.data
    # the previously seen gradient for each sample
    cdef double* gradient_memory = <double*> gradient_memory_init.data

    # the cumulative sums needed for JIT params
    cdef np.ndarray[double, ndim=1] cumulative_sums_array = \
        np.empty(n_samples, dtype=np.double, order="c")
    cdef double* cumulative_sums = <double*> cumulative_sums_array.data

    # the index for the last time this feature was updated
    cdef np.ndarray[int, ndim=1] feature_hist_array = \
        np.zeros(n_features, dtype=np.int32, order="c")
    cdef int* feature_hist = <int*> feature_hist_array.data

    # the previous weights to use to compute stopping criteria
    cdef np.ndarray[double, ndim=2] previous_weights_array = \
        np.zeros((n_features, n_classes), dtype=np.double, order="c")
    cdef double* previous_weights = <double*> previous_weights_array.data

    cdef np.ndarray[double, ndim=1] prediction_array = \
        np.zeros(n_classes, dtype=np.double, order="c")
    cdef double* prediction = <double*> prediction_array.data

    cdef np.ndarray[double, ndim=1] gradient_array = \
        np.zeros(n_classes, dtype=np.double, order="c")
    cdef double* gradient = <double*> gradient_array.data

    # Bias correction term in saga
    cdef double gradient_correction

    # the scalar used for multiplying z
    cdef double wscale = 1.0

    # the cumulative sums for each iteration for the sparse implementation
    cumulative_sums[0] = 0.0

    # the multipliative scale needed for JIT params
    cdef np.ndarray[double, ndim=1] cumulative_sums_prox_array
    cdef double* cumulative_sums_prox

    cdef bint prox = beta > 0 and saga

    # Loss function to optimize
    cdef LossFunction loss
    # Wether the loss function is multinomial
    cdef bint multinomial = False
    # Multinomial loss function
    cdef MultinomialLogLoss multiloss

    if loss_function == "multinomial":
        multinomial = True
        multiloss = MultinomialLogLoss()
    elif loss_function == "log":
        loss = Log()
    elif loss_function == "squared":
        loss = SquaredLoss()
    else:
        raise ValueError("Invalid loss parameter: got %s instead of "
                         "one of ('log', 'squared', 'multinomial')"
                         % loss_function)

    if prox:
        cumulative_sums_prox_array = np.empty(n_samples,
                                              dtype=np.double, order="c")
        cumulative_sums_prox = <double*> cumulative_sums_prox_array.data
    else:
        cumulative_sums_prox = NULL

    with nogil:
        start_time = time(NULL)
        for n_iter in range(max_iter):
            for sample_itr in range(n_samples):
                # extract a random sample
                sample_ind = dataset.random(&x_data_ptr, &x_ind_ptr, &xnnz,
                                              &y, &sample_weight)

                # cached index for gradient_memory
                s_idx = sample_ind * n_classes

                # update the number of samples seen and the seen array
                if seen[sample_ind] == 0:
                    num_seen += 1
                    seen[sample_ind] = 1

                # make the weight updates
                if sample_itr > 0:
                   lagged_update(weights, wscale, xnnz,
                          n_samples, n_classes, sample_itr,
                          cumulative_sums,
                          cumulative_sums_prox,
                          feature_hist,
                          prox,
                          sum_gradient,
                          x_ind_ptr,
                          False,
                          n_iter)

                # find the current prediction
                predict_sample(x_data_ptr, x_ind_ptr, xnnz, weights, wscale,
                               intercept, prediction, n_classes)

                # compute the gradient for this sample, given the prediction
                if multinomial:
                    multiloss._dloss(prediction, y, n_classes, sample_weight,
                                     gradient)
                else:
                    gradient[0] = loss._dloss(prediction[0], y) * sample_weight

                # L2 regularization by simply rescaling the weights
                wscale *= wscale_update

                # make the updates to the sum of gradients
                for j in range(xnnz):
                    feature_ind = x_ind_ptr[j]
                    val = x_data_ptr[j]
                    f_idx = feature_ind * n_classes
                    for class_ind in range(n_classes):
                        gradient_correction = \
                            val * (gradient[class_ind] -
                                   gradient_memory[s_idx + class_ind])
                        if saga:
                            weights[f_idx + class_ind] -= \
                                (gradient_correction * step_size
                                 * (1 - 1. / num_seen) / wscale)
                        sum_gradient[f_idx + class_ind] += gradient_correction

                # fit the intercept
                if fit_intercept:
                    for class_ind in range(n_classes):
                        gradient_correction = (gradient[class_ind] -
                                               gradient_memory[s_idx + class_ind])
                        intercept_sum_gradient[class_ind] += gradient_correction
                        gradient_correction *= step_size * (1. - 1. / num_seen)
                        if saga:
                            intercept[class_ind] -= \
                                (step_size * intercept_sum_gradient[class_ind] /
                                 num_seen * intercept_decay) + gradient_correction
                        else:
                            intercept[class_ind] -= \
                                (step_size * intercept_sum_gradient[class_ind] /
                                 num_seen * intercept_decay)

                        # check to see that the intercept is not inf or NaN
                        if not skl_isfinite(intercept[class_ind]):
                            with gil:
                                raise_infinite_error(n_iter)

                # update the gradient memory for this sample
                for class_ind in range(n_classes):
                    gradient_memory[s_idx + class_ind] = gradient[class_ind]

                if sample_itr == 0:
                    cumulative_sums[0] = step_size / (wscale * num_seen)
                    if prox:
                        cumulative_sums_prox[0] = step_size * beta / wscale
                else:
                    cumulative_sums[sample_itr] = \
                        (cumulative_sums[sample_itr - 1] +
                         step_size / (wscale * num_seen))
                    if prox:
                        cumulative_sums_prox[sample_itr] = \
                        (cumulative_sums_prox[sample_itr - 1] +
                             step_size * beta / wscale)
                # If wscale gets too small, we need to reset the scale.
                if wscale < 1e-9:
                    if verbose:
                        with gil:
                            print("rescaling...")
                    wscale = scale_weights(
                        weights, wscale, n_features, n_samples, n_classes,
                        sample_itr, cumulative_sums,
                        cumulative_sums_prox,
                        feature_hist,
                        prox, sum_gradient, n_iter)

            # we scale the weights every n_samples iterations and reset the
            # just-in-time update system for numerical stability.
            wscale = scale_weights(weights, wscale, n_features, n_samples,
                                   n_classes, n_samples - 1, cumulative_sums,
                                   cumulative_sums_prox,
                                   feature_hist,
                                   prox, sum_gradient, n_iter)

            # check if the stopping criteria is reached
            max_change = 0.0
            max_weight = 0.0
            for idx in range(n_features * n_classes):
                max_weight = fmax(max_weight, fabs(weights[idx]))
                max_change = fmax(max_change,
                                  fabs(weights[idx] -
                                       previous_weights[idx]))
                previous_weights[idx] = weights[idx]
            if ((max_weight != 0 and max_change / max_weight <= tol)
                or max_weight == 0 and max_change == 0):
                if verbose:
                    end_time = time(NULL)
                    with gil:
                        print("convergence after %d epochs took %d seconds" %
                              (n_iter + 1, end_time - start_time))
                break
            elif verbose:
                printf('Epoch %d, change: %.8f\n', n_iter + 1,
                                                  max_change / max_weight)
    n_iter += 1

    if verbose and n_iter >= max_iter:
        end_time = time(NULL)
        print(("max_iter reached after %d seconds") %
              (end_time - start_time))

    return num_seen, n_iter


cdef void raise_infinite_error(int n_iter):
    raise ValueError("Floating-point under-/overflow occurred at "
                     "epoch #%d. Lowering the step_size or "
                     "scaling the input data with StandardScaler "
                     "or MinMaxScaler might help." % (n_iter + 1))


cdef double scale_weights(double* weights, double wscale, int n_features,
                          int n_samples, int n_classes, int sample_itr,
                          double* cumulative_sums,
                          double* cumulative_sums_prox,
                          int* feature_hist,
                          bint prox,
                          double* sum_gradient,
                          int n_iter) nogil:
    """Scale the weights with wscale for numerical stability.

    wscale = (1 - step_size * alpha) ** (n_iter * n_samples + sample_itr)
    can become very small, so we reset it every n_samples iterations to 1.0 for
    numerical stability. To be able to scale, we first need to update every
    coefficients and reset the just-in-time update system.
    This also limits the size of `cumulative_sums`.
    """

    lagged_update(weights, wscale, n_features,
                          n_samples, n_classes, sample_itr + 1,
                          cumulative_sums,
                          cumulative_sums_prox,
                          feature_hist,
                          prox,
                          sum_gradient,
                          NULL,
                          True,
                          n_iter)
    # reset wscale to 1.0
    return 1.0


cdef void lagged_update(double* weights, double wscale, int xnnz,
                          int n_samples, int n_classes, int sample_itr,
                          double* cumulative_sums,
                          double* cumulative_sums_prox,
                          int* feature_hist,
                          bint prox,
                          double* sum_gradient,
                          int* x_ind_ptr,
                          bint reset,
                          int n_iter) nogil:
    """Hard perform the JIT updates for non-zero features of present sample.
     
    The updates that awaits are kept in memory using cumulative_sums,
    cumulative_sums_prox, wscale and feature_hist. See original SAGA paper
    (Defazio et al. 2014) for details. If reset=True, we also reset wscale to
    1 (this is done at the end of each epoch).
    """
    cdef int feature_ind, class_ind, idx, f_idx, lagged_ind, last_update_ind
    cdef double cum_sum, grad_step, prox_step
    for feature_ind in range(xnnz):
        if not reset:
            feature_ind = x_ind_ptr[feature_ind]
        f_idx = feature_ind * n_classes

        cum_sum = cumulative_sums[sample_itr - 1]
        if prox:
            cum_sum_prox = cumulative_sums_prox[sample_itr - 1]
        if feature_hist[feature_ind] != 0:
            cum_sum -= cumulative_sums[feature_hist[feature_ind] - 1]
            if prox:
                cum_sum_prox -= cumulative_sums_prox[feature_hist[feature_ind] - 1]
        if not prox:
            for class_ind in range(n_classes):
                idx = f_idx + class_ind
                weights[idx] -= cum_sum * sum_gradient[idx]
                if reset:
                    weights[idx] *= wscale
                    if not skl_isfinite(weights[idx]):
                        with gil:
                            raise_infinite_error(n_iter)
        else:
            for class_ind in range(n_classes):
                idx = f_idx + class_ind
                if fabs(sum_gradient[idx] * cum_sum) < cum_sum_prox:
                    # In this case, we can perform all the gradient steps and
                    # all the proximal steps in this order, which is more
                    # efficient than unrolling all the lagged updates.
                    # Idea taken from scikit-learn-contrib/lightning.
                    weights[idx] -= cum_sum * sum_gradient[idx]
                    weights[idx] = _soft_thresholding(weights[idx],
                                                      cum_sum_prox)
                else:
                    last_update_ind = feature_hist[feature_ind] - 1
                    if last_update_ind == -1:
                        last_update_ind = sample_itr - 1
                    for lagged_ind in range(sample_itr - 1,
                                   last_update_ind - 1, -1):
                        if lagged_ind > 0:
                            grad_step = (cumulative_sums[lagged_ind]
                               - cumulative_sums[lagged_ind - 1])
                            prox_step = (cumulative_sums_prox[lagged_ind]
                               - cumulative_sums_prox[lagged_ind - 1])
                        else:
                            grad_step = cumulative_sums[lagged_ind]
                            prox_step = cumulative_sums_prox[lagged_ind]
                        weights[idx] -= sum_gradient[idx] * grad_step
                        weights[idx] = _soft_thresholding(weights[idx],
                                                          prox_step)

                if reset:
                    weights[idx] *= wscale
                    # check to see that the weight is not inf or NaN
                    if not skl_isfinite(weights[idx]):
                        with gil:
                            raise_infinite_error(n_iter)
        if reset:
            feature_hist[feature_ind] = sample_itr % n_samples
        else:
            feature_hist[feature_ind] = sample_itr

    if reset:
        cumulative_sums[sample_itr - 1] = 0.0
        if prox:
            cumulative_sums_prox[sample_itr - 1] = 0.0


cdef void predict_sample(double* x_data_ptr, int* x_ind_ptr, int xnnz,
                         double* w_data_ptr, double wscale, double* intercept,
                         double* prediction, int n_classes) nogil:
    """Compute the prediction given sparse sample x and dense weight w.

    Parameters
    ----------
    x_data_ptr : pointer
        Pointer to the data of the sample x

    x_ind_ptr : pointer
        Pointer to the indices of the sample  x

    xnnz : int
        Number of non-zero element in the sample  x

    w_data_ptr : pointer
        Pointer to the data of the weights w

    wscale : double
        Scale of the weights w

    intercept : pointer
        Pointer to the intercept

    prediction : pointer
        Pointer to store the resulting prediction

    n_classes : int
        Number of classes in multinomial case. Equals 1 in binary case.

    """
    cdef int feature_ind, class_ind, j
    cdef double innerprod

    for class_ind in range(n_classes):
        innerprod = 0.0
        # Compute the dot product only on non-zero elements of x
        for j in range(xnnz):
            feature_ind = x_ind_ptr[j]
            innerprod += (w_data_ptr[feature_ind * n_classes + class_ind] *
                          x_data_ptr[j])

        prediction[class_ind] = wscale * innerprod + intercept[class_ind]
