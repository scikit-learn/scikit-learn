# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3

# Design:
# See https://github.com/scikit-learn/scikit-learn/issues/15123 for reasons.
# a) Merge link functions into loss functions for speed and numerical
#    stability, i.e. use raw_prediction instead of y_pred in signature.
# b) Pure C functions (nogil) calculate single points (single sample)
# c) Wrap C functions in a loop to get Python functions operating on ndarrays.
#   - Write loops manually.
#     Reason: There is still some performance overhead when using a wrapper
#     function "wrap" that carries out the loop and gets as argument a function
#     pointer to one of the C functions from b), e.g.
#     wrap(closs_half_poisson, y_true, ...)
#   - Pass n_threads as argument to prange and propagate option to all callers.
# d) Provide classes (Cython extension types) per loss in order to have
#    semantical structured objects.
#    - Member function for single points just call the C function from b).
#      These are used e.g. in SGD `_plain_sgd`.
#    - Member functions operating on ndarrays looping, see c), over calls to C
#      functions from b).
# e) Provide convenience Python classes that inherit from these extension types
#    elsewhere (see loss.py)
#    - Example: loss.gradient calls extension_type._gradient but does some
#      input checking like None -> np.empty().
#
# Note: We require 1-dim ndarrays to be contiguous.
# TODO: Use const memoryviews with fused types with Cython 3.0 where
#       appropriate (arguments marked by "# IN")

cimport cython
from cython.parallel import parallel, prange
import numpy as np
cimport numpy as np

from libc.math cimport exp, fabs, log, log1p
from libc.stdlib cimport malloc, free

np.import_array()


# -------------------------------------
# Helper functions
# -------------------------------------
# Numerically stable version of log(1 + exp(x)) for double precision
# See https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
cdef inline double log1pexp(double x) nogil:
    if x <= -37:
        return exp(x)
    elif x <= 18:
        return log1p(exp(x))
    elif x <= 33.3:
        return x + exp(-x)
    else:
        return x


cdef inline void sum_exp_minus_max(
    const int i,
    Y_DTYPE_C[:, :] raw_prediction,  # IN
    Y_DTYPE_C *p                     # OUT
) nogil:
    # Store p[k] = exp(raw_prediction_i_k - max_value) for k = 0 to n_classes-1
    #       p[-2] = max(raw_prediction_i_k, k = 0 to n_classes-1)
    #       p[-1] = sum(p[k], k = 0 to n_classes-1) = sum of exponentials
    # len(p) must be n_classes + 2
    # Notes:
    # - Using "by reference" arguments doesn't work well, therefore we use a
    #   longer p, see https://github.com/cython/cython/issues/1863
    # - i needs to be passed (and stays constant) because otherwise Cython does
    #   not generate optimal code, see
    #   https://github.com/scikit-learn/scikit-learn/issues/17299
    # - We do not normalize p by calculating p[k] = p[k] / sum_exps.
    #   This helps to save one loop over k.
    cdef:
        int k
        int n_classes = raw_prediction.shape[1]
        double max_value = raw_prediction[i, 0]
        double sum_exps = 0
    for k in range(1, n_classes):
        # Compute max value of array for numerical stability
        if max_value < raw_prediction[i, k]:
            max_value = raw_prediction[i, k]

    for k in range(n_classes):
        p[k] = exp(raw_prediction[i, k] - max_value)
        sum_exps += p[k]

    p[n_classes] = max_value     # same as p[-2]
    p[n_classes + 1] = sum_exps  # same as p[-1]


# -------------------------------------
# Single point inline C functions
# -------------------------------------
# Half Squared Error
cdef inline double closs_half_squared_error(
    double y_true,
    double raw_prediction
) nogil:
    return 0.5 * (raw_prediction - y_true) * (raw_prediction - y_true)


cdef inline double cgradient_half_squared_error(
    double y_true,
    double raw_prediction
) nogil:
    return raw_prediction - y_true


cdef inline double2 cgrad_hess_half_squared_error(
    double y_true,
    double raw_prediction
) nogil:
    cdef double2 gh
    gh.val1 = raw_prediction - y_true  # gradient
    gh.val2 = 1.                       # hessian
    return gh


# Absolute Error
cdef inline double closs_absolute_error(
    double y_true,
    double raw_prediction
) nogil:
    return fabs(raw_prediction - y_true)


cdef inline double cgradient_absolute_error(
    double y_true,
    double raw_prediction
) nogil:
    return 1. if raw_prediction > y_true else -1.


cdef inline double2 cgrad_hess_absolute_error(
    double y_true,
    double raw_prediction
) nogil:
    cdef double2 gh
    # Note that exact hessian = 0 almost everywhere. Optimization routines like
    # in HGBT, however, need a hessian > 0. Therefore, we assign 1.
    gh.val1 = 1. if raw_prediction > y_true else -1.  # gradient
    gh.val2 = 1.                                      # hessian
    return gh


# Quantile Loss / Pinball Loss
cdef inline double closs_pinball_loss(
    double y_true,
    double raw_prediction,
    double quantile
) nogil:
    return (quantile * (y_true - raw_prediction) if y_true >= raw_prediction
            else (1. - quantile) * (raw_prediction - y_true))


cdef inline double cgradient_pinball_loss(
    double y_true,
    double raw_prediction,
    double quantile
) nogil:
    return -quantile if y_true >=raw_prediction else 1. - quantile


cdef inline double2 cgrad_hess_pinball_loss(
    double y_true,
    double raw_prediction,
    double quantile
) nogil:
    cdef double2 gh
    # Note that exact hessian = 0 almost everywhere. Optimization routines like
    # in HGBT, however, need a hessian > 0. Therefore, we assign 1.
    gh.val1 = -quantile if y_true >=raw_prediction else 1. - quantile  # gradient
    gh.val2 = 1.                                                       # hessian
    return gh


# Half Poisson Deviance with Log-Link, dropping constant terms
cdef inline double closs_half_poisson(
    double y_true,
    double raw_prediction
) nogil:
    return exp(raw_prediction) - y_true * raw_prediction


cdef inline double cgradient_half_poisson(
    double y_true,
    double raw_prediction
) nogil:
    # y_pred - y_true
    return exp(raw_prediction) - y_true


cdef inline double2 closs_grad_half_poisson(
    double y_true,
    double raw_prediction
) nogil:
    cdef double2 lg
    lg.val2 = exp(raw_prediction)                # used as temporary
    lg.val1 = lg.val2 - y_true * raw_prediction  # loss
    lg.val2 -= y_true                            # gradient
    return lg


cdef inline double2 cgrad_hess_half_poisson(
    double y_true,
    double raw_prediction
) nogil:
    cdef double2 gh
    gh.val2 = exp(raw_prediction)  # hessian
    gh.val1 = gh.val2 - y_true     # gradient
    return gh


# Half Gamma Deviance with Log-Link, dropping constant terms
cdef inline double closs_half_gamma(
    double y_true,
    double raw_prediction
) nogil:
    return raw_prediction + y_true * exp(-raw_prediction)


cdef inline double cgradient_half_gamma(
    double y_true,
    double raw_prediction
) nogil:
    return 1. - y_true * exp(-raw_prediction)


cdef inline double2 closs_grad_half_gamma(
    double y_true,
    double raw_prediction
) nogil:
    cdef double2 lg
    lg.val2 = exp(-raw_prediction)               # used as temporary
    lg.val1 = raw_prediction + y_true * lg.val2  # loss
    lg.val2 = 1. - y_true * lg.val2              # gradient
    return lg


cdef inline double2 cgrad_hess_half_gamma(
    double y_true,
    double raw_prediction
) nogil:
    cdef double2 gh
    gh.val2 = exp(-raw_prediction)   # used as temporary
    gh.val1 = 1. - y_true * gh.val2  # gradient
    gh.val2 *= y_true                # hessian
    return gh


# Half Tweedie Deviance with Log-Link, dropping constant terms
# Note that by dropping constants this is no longer smooth in parameter power.
cdef inline double closs_half_tweedie(
    double y_true,
    double raw_prediction,
    double power
) nogil:
    if power == 0.:
        return closs_half_squared_error(y_true, exp(raw_prediction))
    elif power == 1.:
        return closs_half_poisson(y_true, raw_prediction)
    elif power == 2.:
        return closs_half_gamma(y_true, raw_prediction)
    else:
        return (exp((2. - power) * raw_prediction) / (2. - power)
                - y_true * exp((1. - power) * raw_prediction) / (1. - power))


cdef inline double cgradient_half_tweedie(
    double y_true,
    double raw_prediction,
    double power
) nogil:
    cdef double exp1
    if power == 0.:
        exp1 = exp(raw_prediction)
        return exp1 * (exp1 - y_true)
    elif power == 1.:
        return cgradient_half_poisson(y_true, raw_prediction)
    elif power == 2.:
        return cgradient_half_gamma(y_true, raw_prediction)
    else:
        return (exp((2. - power) * raw_prediction)
                - y_true * exp((1. - power) * raw_prediction))


cdef inline double2 closs_grad_half_tweedie(
    double y_true,
    double raw_prediction,
    double power
) nogil:
    cdef double2 lg
    cdef double exp1, exp2
    if power == 0.:
        exp1 = exp(raw_prediction)
        lg.val1 = closs_half_squared_error(y_true, exp1)  # loss
        lg.val2 = exp1 * (exp1 - y_true)                  # gradient
    elif power == 1.:
        return closs_grad_half_poisson(y_true, raw_prediction)
    elif power == 2.:
        return closs_grad_half_gamma(y_true, raw_prediction)
    else:
        exp1 = exp((1. - power) * raw_prediction)
        exp2 = exp((2. - power) * raw_prediction)
        lg.val1 = exp2 / (2. - power) - y_true * exp1 / (1. - power)  # loss
        lg.val2 = exp2 - y_true * exp1                                # gradient
    return lg


cdef inline double2 cgrad_hess_half_tweedie(
    double y_true,
    double raw_prediction,
    double power
) nogil:
    cdef double2 gh
    cdef double exp1, exp2
    if power == 0.:
        exp1 = exp(raw_prediction)
        gh.val1 = exp1 * (exp1 - y_true)      # gradient
        gh.val2 = exp1 * (2 * exp1 - y_true)  # hessian
    elif power == 1.:
        return cgrad_hess_half_poisson(y_true, raw_prediction)
    elif power == 2.:
        return cgrad_hess_half_gamma(y_true, raw_prediction)
    else:
        exp1 = exp((1. - power) * raw_prediction)
        exp2 = exp((2. - power) * raw_prediction)
        gh.val1 = exp2 - y_true * exp1                                # gradient
        gh.val2 = (2. - power) * exp2 - (1. - power) * y_true * exp1  # hessian
    return gh


# Binary cross entropy aka log-loss
cdef inline double closs_binary_crossentropy(
    double y_true,
    double raw_prediction
) nogil:
    # log1p(exp(raw_prediction)) - y_true * raw_prediction
    return log1pexp(raw_prediction) - y_true * raw_prediction


cdef inline double cgradient_binary_crossentropy(
    double y_true,
    double raw_prediction
) nogil:
    # y_pred - y_true = expit(raw_prediction) - y_true
    # Numerically more stable, see
    # http://fa.bianp.net/blog/2019/evaluate_logistic/
    #     if raw_prediction < 0:
    #         exp_tmp = exp(raw_prediction)
    #         return ((1 - y_true) * exp_tmp - y_true) / (1 + exp_tmp)
    #     else:
    #         exp_tmp = exp(-raw_prediction)
    #         return ((1 - y_true) - y_true * exp_tmp) / (1 + exp_tmp)
    # Note that optimal speed would be achieved, at the cost of precision, by
    #     return expit(raw_prediction) - y_true
    # i.e. no if else, and an own inline implemention of expit instead of
    #     from scipy.special.cython_special cimport expit
    # The case distinction raw_prediction < 0 in the stable implementation
    # does not provide significant better precision. Therefore we go without
    # it.
    cdef double exp_tmp
    exp_tmp = exp(-raw_prediction)
    return ((1 - y_true) - y_true * exp_tmp) / (1 + exp_tmp)


cdef inline double2 closs_grad_binary_crossentropy(
    double y_true,
    double raw_prediction
) nogil:
    cdef double2 lg
    if raw_prediction <= 0:
        lg.val2 = exp(raw_prediction)  # used as temporary
        if raw_prediction <= -37:
            lg.val1 = lg.val2 - y_true * raw_prediction              # loss
        else:
            lg.val1 = log1p(lg.val2) - y_true * raw_prediction       # loss
        lg.val2 = ((1 - y_true) * lg.val2 - y_true) / (1 + lg.val2)  # gradient
    else:
        lg.val2 = exp(-raw_prediction)  # used as temporary
        if raw_prediction <= 18:
            # log1p(exp(x)) = log(1 + exp(x)) = x + log1p(exp(-x))
            lg.val1 = log1p(lg.val2) + (1 - y_true) * raw_prediction  # loss
        else:
            lg.val1 = lg.val2 + (1 - y_true) * raw_prediction         # loss
        lg.val2 = ((1 - y_true) - y_true * lg.val2) / (1 + lg.val2)   # gradient
    return lg


cdef inline double2 cgrad_hess_binary_crossentropy(
    double y_true,
    double raw_prediction
) nogil:
    # with y_pred = expit(raw)
    # hessian = y_pred * (1 - y_pred) = exp(raw) / (1 + exp(raw))**2
    #                                 = exp(-raw) / (1 + exp(-raw))**2
    cdef double2 gh
    gh.val2 = exp(-raw_prediction)  # used as temporary
    gh.val1 = ((1 - y_true) - y_true * gh.val2) / (1 + gh.val2)  # gradient
    gh.val2 = gh.val2 / (1 + gh.val2)**2                         # hessian
    return gh


# ---------------------------------------------------
# Extension Types for Loss Functions of 1-dim targets
# ---------------------------------------------------
cdef class cLossFunction:
    """Base class for convex loss functions."""

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        """Compute the loss for a single sample.

        Parameters
        ----------
        y_true : double
            Observed, true target value.
        raw_prediction : double
            Raw prediction value (in link space).

        Returns
        -------
        double
            The loss evaluated at `y_true` and `raw_prediction`.
        """
        pass

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        """Compute gradient of loss w.r.t. raw_prediction for a single sample.

        Parameters
        ----------
        y_true : double
            Observed, true target value.
        raw_prediction : double
            Raw prediction value (in link space).

        Returns
        -------
        double
            The derivative of the loss function w.r.t. `raw_prediction`.
        """
        pass

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        """Compute gradient and hessian.

        Gradient and hessian of loss w.r.t. raw_prediction for a single sample.

        This is usually diagonal in raw_prediction_i and raw_prediction_j.
        Therefore, we return the diagonal element i=j.

        For a loss with a non-canonical link, this might implement the diagonal
        of the Fisher matrix (=expected hessian) instead of the hessian.

        Parameters
        ----------
        y_true : double
            Observed, true target value.
        raw_prediction : double
            Raw prediction value (in link space).

        Returns
        -------
        grad_hess_pair
            Gradient and hessian of the loss function w.r.t. `raw_prediction`.
        """
        pass

    # Note: With Cython 3.0, fused types can be used together with const:
    #       const Y_DTYPE_C double[::1] y_true
    # See release notes 3.0.0 alpha1
    # https://cython.readthedocs.io/en/latest/src/changes.html#alpha-1-2020-04-12
    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,          # IN
        Y_DTYPE_C[::1] raw_prediction,  # IN
        Y_DTYPE_C[::1] sample_weight,   # IN
        G_DTYPE_C[::1] loss,            # OUT
        int n_threads=1
    ):
        """Compute the pointwise loss value for each input.

        Parameters
        ----------
        y_true : array of shape (n_samples,)
            Observed, true target values.
        raw_prediction : array of shape (n_samples,)
            Raw prediction values (in link space).
        sample_weight : array of shape (n_samples,) or None
            Sample weights.
        loss : array of shape (n_samples,)
            A location into which the result is stored.
        n_threads : int
            Might use openmp thread parallelism.

        Returns
        -------
        loss : array of shape (n_samples,)
            Element-wise loss function.
        """
        pass

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,          # IN
        Y_DTYPE_C[::1] raw_prediction,  # IN
        Y_DTYPE_C[::1] sample_weight,   # IN
        G_DTYPE_C[::1] gradient,        # OUT
        int n_threads=1
    ):
        """Compute gradient of loss w.r.t raw_prediction for each input.

        Parameters
        ----------
        y_true : array of shape (n_samples,)
            Observed, true target values.
        raw_prediction : array of shape (n_samples,)
            Raw prediction values (in link space).
        sample_weight : array of shape (n_samples,) or None
            Sample weights.
        gradient : array of shape (n_samples,)
            A location into which the result is stored.
        n_threads : int
            Might use openmp thread parallelism.

        Returns
        -------
        gradient : array of shape (n_samples,)
            Element-wise gradients.
        """
        pass

    def _loss_gradient(
        self,
        Y_DTYPE_C[::1] y_true,          # IN
        Y_DTYPE_C[::1] raw_prediction,  # IN
        Y_DTYPE_C[::1] sample_weight,   # IN
        G_DTYPE_C[::1] loss,            # OUT
        G_DTYPE_C[::1] gradient,        # OUT
        int n_threads=1
    ):
        """Compute loss and gradient of loss w.r.t raw_prediction.

        Parameters
        ----------
        y_true : array of shape (n_samples,)
            Observed, true target values.
        raw_prediction : array of shape (n_samples,)
            Raw prediction values (in link space).
        sample_weight : array of shape (n_samples,) or None
            Sample weights.
        loss : array of shape (n_samples,) or None
            A location into which the element-wise loss is stored.
        gradient : array of shape (n_samples,)
            A location into which the gradient is stored.
        n_threads : int
            Might use openmp thread parallelism.

        Returns
        -------
        loss : array of shape (n_samples,)
            Element-wise loss function.

        gradient : array of shape (n_samples,)
            Element-wise gradients.
        """
        self._loss(y_true, raw_prediction, sample_weight, loss,
                            n_threads)
        self._gradient(y_true, raw_prediction, sample_weight, gradient,
                      n_threads)
        return np.asarray(loss), np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,          # IN
        Y_DTYPE_C[::1] raw_prediction,  # IN
        Y_DTYPE_C[::1] sample_weight,   # IN
        G_DTYPE_C[::1] gradient,        # OUT
        G_DTYPE_C[::1] hessian,         # OUT
        int n_threads=1
    ):
        """Compute gradient and hessian of loss w.r.t raw_prediction.

        Parameters
        ----------
        y_true : array of shape (n_samples,)
            Observed, true target values.
        raw_prediction : array of shape (n_samples,)
            Raw prediction values (in link space).
        sample_weight : array of shape (n_samples,) or None
            Sample weights.
        gradient : array of shape (n_samples,)
            A location into which the gradient is stored.
        hessian : array of shape (n_samples,)
            A location into which the hessian is stored.
        n_threads : int
            Might use openmp thread parallelism.

        Returns
        -------
        gradient : array of shape (n_samples,)
            Element-wise gradients.

        hessian : array of shape (n_samples,)
            Element-wise hessians.
        """
        pass


cdef class cHalfSquaredError(cLossFunction):
    """Half Squared Error with identity link.

    Domain:
    y_true and y_pred all real numbers

    Link:
    y_pred = raw_prediction
    """

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        return closs_half_squared_error(y_true, raw_prediction)

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        return cgradient_half_squared_error(y_true, raw_prediction)

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        return cgrad_hess_half_squared_error(y_true, raw_prediction)

    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = closs_half_squared_error(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = (
                    sample_weight[i]
                    * closs_half_squared_error(y_true[i], raw_prediction[i])
                )

        return np.asarray(loss)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = cgradient_half_squared_error(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = (
                    sample_weight[i]
                    * cgradient_half_squared_error(y_true[i], raw_prediction[i])
                )

        return np.asarray(gradient)


    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        G_DTYPE_C[::1] hessian,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_squared_error(y_true[i], raw_prediction[i])
                gradient[i] = dbl2.val1
                hessian[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_squared_error(y_true[i], raw_prediction[i])
                gradient[i] = sample_weight[i] * dbl2.val1
                hessian[i] = sample_weight[i] * dbl2.val2

        return np.asarray(gradient), np.asarray(hessian)


cdef class cAbsoluteError(cLossFunction):
    """Absolute Error with identity link.

    Domain:
    y_true and y_pred all real numbers

    Link:
    y_pred = raw_prediction
    """

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        return closs_absolute_error(y_true, raw_prediction)

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        return cgradient_absolute_error(y_true, raw_prediction)

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        return cgrad_hess_absolute_error(y_true, raw_prediction)

    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = closs_absolute_error(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = (sample_weight[i]
                    * closs_absolute_error(y_true[i], raw_prediction[i]))

        return np.asarray(loss)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = cgradient_absolute_error(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = (
                    sample_weight[i]
                    * cgradient_absolute_error(y_true[i], raw_prediction[i])
                )

        return np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        G_DTYPE_C[::1] hessian,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_absolute_error(y_true[i], raw_prediction[i])
                gradient[i] = dbl2.val1
                hessian[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_absolute_error(y_true[i], raw_prediction[i])
                gradient[i] = sample_weight[i] * dbl2.val1
                hessian[i] = sample_weight[i] * dbl2.val2

        return np.asarray(gradient), np.asarray(hessian)


cdef class cPinballLoss(cLossFunction):
    """Quantile Loss aka Pinball Loss with identity link.

    Domain:
    y_true and y_pred all real numbers
    quantile in (0, 1)

    Link:
    y_pred = raw_prediction

    Note: 2 * cPinballLoss(quantile=0.5) equals cAbsoluteError()
    """

    def __init__(self, quantile):
        self.quantile = quantile

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        return closs_pinball_loss(y_true, raw_prediction, self.quantile)

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        return cgradient_pinball_loss(y_true, raw_prediction, self.quantile)

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        return cgrad_hess_pinball_loss(y_true, raw_prediction, self.quantile)

    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = closs_pinball_loss(y_true[i], raw_prediction[i], self.quantile)
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = (
                    sample_weight[i]
                    * closs_pinball_loss(y_true[i], raw_prediction[i], self.quantile)
                )

        return np.asarray(loss)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = cgradient_pinball_loss(
                    y_true[i], raw_prediction[i], self.quantile
                )
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = (
                    sample_weight[i]
                    * cgradient_pinball_loss(y_true[i], raw_prediction[i], self.quantile)
                )

        return np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        G_DTYPE_C[::1] hessian,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_pinball_loss(
                    y_true[i], raw_prediction[i], self.quantile
                )
                gradient[i] = dbl2.val1
                hessian[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_pinball_loss(
                    y_true[i], raw_prediction[i], self.quantile
                )
                gradient[i] = sample_weight[i] * dbl2.val1
                hessian[i] = sample_weight[i] * dbl2.val2

        return np.asarray(gradient), np.asarray(hessian)


cdef class cHalfPoissonLoss(cLossFunction):
    """Half Poisson deviance loss with log-link.

    Domain:
    y_true in non-negative real numbers
    y_pred in positive real numbers

    Link:
    y_pred = exp(raw_prediction)

    Half Poisson deviance with log-link is
        y_true * log(y_true/y_pred) + y_pred - y_true
        = y_true * log(y_true) - y_true * raw_prediction
          + exp(raw_prediction) - y_true

    Dropping constant terms, this gives:
        exp(raw_prediction) - y_true * raw_prediction
    """

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        return closs_half_poisson(y_true, raw_prediction)

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        return cgradient_half_poisson(y_true, raw_prediction)

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        return cgrad_hess_half_poisson(y_true, raw_prediction)

    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = closs_half_poisson(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = (
                    sample_weight[i]
                    * closs_half_poisson(y_true[i], raw_prediction[i])
                )

        return np.asarray(loss)

    def _loss_gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_half_poisson(y_true[i], raw_prediction[i])
                loss[i] = dbl2.val1
                gradient[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_half_poisson(y_true[i], raw_prediction[i])
                loss[i] = sample_weight[i] * dbl2.val1
                gradient[i] = sample_weight[i] * dbl2.val2

        return np.asarray(loss), np.asarray(gradient)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = cgradient_half_poisson(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = (
                    sample_weight[i]
                    * cgradient_half_poisson(y_true[i], raw_prediction[i])
                )

        return np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        G_DTYPE_C[::1] hessian,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_poisson(y_true[i], raw_prediction[i])
                gradient[i] = dbl2.val1
                hessian[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_poisson(y_true[i], raw_prediction[i])
                gradient[i] = sample_weight[i] * dbl2.val1
                hessian[i] = sample_weight[i] * dbl2.val2

        return np.asarray(gradient), np.asarray(hessian)


cdef class cHalfGammaLoss(cLossFunction):
    """Half Gamma deviance loss with log-link.

    Domain:
    y_true and y_pred in positive real numbers

    Link:
    y_pred = exp(raw_prediction)

    Half Gamma deviance with log-link is
        log(y_pred/y_true) + y_true/y_pred - 1
        = raw_prediction - log(y_true) + y_true * exp(-raw_prediction) - 1

    Dropping constant terms, this gives:
        raw_prediction + y_true * exp(-raw_prediction)
    """

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        return closs_half_gamma(y_true, raw_prediction)

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        return cgradient_half_gamma(y_true, raw_prediction)

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        return cgrad_hess_half_gamma(y_true, raw_prediction)

    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = closs_half_gamma(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = (
                    sample_weight[i]
                    * closs_half_gamma(y_true[i], raw_prediction[i])
                )

        return np.asarray(loss)

    def _loss_gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_half_gamma(y_true[i], raw_prediction[i])
                loss[i] = dbl2.val1
                gradient[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_half_gamma(y_true[i], raw_prediction[i])
                loss[i] = sample_weight[i] * dbl2.val1
                gradient[i] = sample_weight[i] * dbl2.val2

        return np.asarray(loss), np.asarray(gradient)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = cgradient_half_gamma(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = (
                    sample_weight[i]
                    * cgradient_half_gamma(y_true[i], raw_prediction[i])
                )

        return np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        G_DTYPE_C[::1] hessian,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_gamma(y_true[i], raw_prediction[i])
                gradient[i] = dbl2.val1
                hessian[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_gamma(y_true[i], raw_prediction[i])
                gradient[i] = sample_weight[i] * dbl2.val1
                hessian[i] = sample_weight[i] * dbl2.val2

        return np.asarray(gradient), np.asarray(hessian)


cdef class cHalfTweedieLoss(cLossFunction):
    """Half Tweedie deviance loss with log-link.

    Domain:
    y_true in real numbers if p <= 0
    y_true in non-negative real numbers if 0 < p < 2
    y_true in positive real numbers if p >= 2
    y_pred and power in positive real numbers

    Link:
    y_pred = exp(raw_prediction)

    Half Tweedie deviance with log-link and p=power is
        max(y_true, 0)**(2-p) / (1-p) / (2-p)
        - y_true * y_pred**(1-p) / (1-p)
        + y_pred**(2-p) / (2-p)
        = max(y_true, 0)**(2-p) / (1-p) / (2-p)
        - y_true * exp((1-p) * raw_prediction) / (1-p)
        + exp((2-p) * raw_prediction) / (2-p)

    Dropping constant terms, this gives:
        exp((2-p) * raw_prediction) / (2-p)
        - y_true * exp((1-p) * raw_prediction) / (1-p)

    Notes:
    - Poisson with p=1 and and Gamma with p=2 have different terms dropped such
      that cHalfTweedieLoss is not continuous in p=power at p=1 and p=2.
    - While the Tweedie distribution only exists for p<=0 or p>=1, the range
      0<p<1 still gives a strictly consistent scoring function for the
      expectation.
    """

    def __init__(self, power):
        self.power = power

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        return closs_half_tweedie(y_true, raw_prediction, self.power)

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        return cgradient_half_tweedie(y_true, raw_prediction, self.power)

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        return cgrad_hess_half_tweedie(y_true, raw_prediction, self.power)

    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = closs_half_tweedie(y_true[i], raw_prediction[i], self.power)
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = (
                    sample_weight[i]
                    * closs_half_tweedie(y_true[i], raw_prediction[i], self.power)
                )

        return np.asarray(loss)

    def _loss_gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_half_tweedie(y_true[i], raw_prediction[i], self.power)
                loss[i] = dbl2.val1
                gradient[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_half_tweedie(y_true[i], raw_prediction[i], self.power)
                loss[i] = sample_weight[i] * dbl2.val1
                gradient[i] = sample_weight[i] * dbl2.val2

        return np.asarray(loss), np.asarray(gradient)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = cgradient_half_tweedie(
                    y_true[i], raw_prediction[i], self.power
                )
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = (
                    sample_weight[i]
                    * cgradient_half_tweedie(y_true[i], raw_prediction[i], self.power)
                )

        return np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        G_DTYPE_C[::1] hessian,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_tweedie(y_true[i], raw_prediction[i], self.power)
                gradient[i] = dbl2.val1
                hessian[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_half_tweedie(y_true[i], raw_prediction[i], self.power)
                gradient[i] = sample_weight[i] * dbl2.val1
                hessian[i] = sample_weight[i] * dbl2.val2

        return np.asarray(gradient), np.asarray(hessian)


cdef class cBinaryCrossEntropy(cLossFunction):
    """BinaryCrossEntropy with logit link.

    Domain:
    y_true in [0, 1]
    y_pred in (0, 1), i.e. boundaries excluded

    Link:
    y_pred = expit(raw_prediction)
    """

    cdef double closs(self, double y_true, double raw_prediction) nogil:
        return closs_binary_crossentropy(y_true, raw_prediction)

    cdef double cgradient(self, double y_true, double raw_prediction) nogil:
        return cgradient_binary_crossentropy(y_true, raw_prediction)

    cdef double2 cgrad_hess(self, double y_true, double raw_prediction) nogil:
        return cgrad_hess_binary_crossentropy(y_true, raw_prediction)

    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = closs_binary_crossentropy(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                loss[i] = (
                    sample_weight[i]
                    * closs_binary_crossentropy(y_true[i], raw_prediction[i])
                )

        return np.asarray(loss)

    def _loss_gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_binary_crossentropy(y_true[i], raw_prediction[i])
                loss[i] = dbl2.val1
                gradient[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = closs_grad_binary_crossentropy(y_true[i], raw_prediction[i])
                loss[i] = sample_weight[i] * dbl2.val1
                gradient[i] = sample_weight[i] * dbl2.val2

        return np.asarray(loss), np.asarray(gradient)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = cgradient_binary_crossentropy(y_true[i], raw_prediction[i])
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                gradient[i] = (
                    sample_weight[i]
                    * cgradient_binary_crossentropy(y_true[i], raw_prediction[i])
                )

        return np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[::1] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] gradient,
        G_DTYPE_C[::1] hessian,
        int n_threads=1
    ):
        cdef:
            int i
            int n_samples = y_true.shape[0]
            double2 dbl2

        if sample_weight is None:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_binary_crossentropy(y_true[i], raw_prediction[i])
                gradient[i] = dbl2.val1
                hessian[i] = dbl2.val2
        else:
            for i in prange(
                n_samples, schedule='static', nogil=True, num_threads=n_threads
            ):
                dbl2 = cgrad_hess_binary_crossentropy(y_true[i], raw_prediction[i])
                gradient[i] = sample_weight[i] * dbl2.val1
                hessian[i] = sample_weight[i] * dbl2.val2

        return np.asarray(gradient), np.asarray(hessian)


cdef class cCategoricalCrossEntropy(cLossFunction):
    """CategoricalCrossEntropy with multinomial logit link.

    Domain:
    y_true in {0, 1, 2, 3, .., n_classes - 1}
    y_pred in (0, 1)**n_classes, i.e. interval with boundaries excluded

    Link:
    y_pred = softmax(raw_prediction)

    Note: Label encoding is built-in, i.e. {0, 1, 2, 3, .., n_classes - 1} is
    mapped to (y_true == k) for k = 0 .. n_classes - 1 which is either 0 or 1.
    """

    # Note that we do not assume memory alignement/contiguity of 2d arrays.
    # There seems to be little benefit in doing so. Benchmarks proofing the
    # opposite are welcome.
    def _loss(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[:, :] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        int n_threads=1
    ):
        cdef:
            int i, k
            int n_samples = y_true.shape[0]
            int n_classes = raw_prediction.shape[1]
            Y_DTYPE_C max_value, sum_exps
            Y_DTYPE_C*  p  # temporary buffer

        # We assume n_samples > n_classes. In this case having the inner loop
        # over n_classes is a good default.
        # TODO: If every memoryview is contiguous and raw_preduction is
        #       f-contiguous, can we write a better algo (loops) to improve
        #       performance?
        if sample_weight is None:
            # inner loop over n_classes
            with nogil, parallel(num_threads=n_threads):
                # Define private buffer variables as each thread might use its
                # own.
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    max_value = p[n_classes]     # p[-2]
                    sum_exps = p[n_classes + 1]  # p[-1]
                    loss[i] = log(sum_exps) + max_value

                    for k in range(n_classes):
                        # label decode y_true
                        if y_true[i] == k:
                            loss[i] -= raw_prediction[i, k]

                free(p)
        else:
            with nogil, parallel(num_threads=n_threads):
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    max_value = raw_prediction[i, 0]
                    max_value = p[n_classes]     # p[-2]
                    sum_exps = p[n_classes + 1]  # p[-1]
                    loss[i] = log(sum_exps) + max_value

                    for k in range(n_classes):
                        # label decode y_true
                        if y_true[i] == k:
                            loss[i] -= raw_prediction[i, k]

                    loss[i] *= sample_weight[i]

                free(p)

        return np.asarray(loss)

    def _loss_gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[:, :] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[::1] loss,
        G_DTYPE_C[:, :] gradient,
        int n_threads=1
    ):
        cdef:
            int i, k
            int n_samples = y_true.shape[0]
            int n_classes = raw_prediction.shape[1]
            Y_DTYPE_C max_value, sum_exps
            Y_DTYPE_C*  p  # temporary buffer

        if sample_weight is None:
            # inner loop over n_classes
            with nogil, parallel(num_threads=n_threads):
                # Define private buffer variables as each thread might use its
                # own.
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    max_value = p[n_classes]  # p[-2]
                    sum_exps = p[n_classes + 1]  # p[-1]
                    loss[i] = log(sum_exps) + max_value

                    for k in range(n_classes):
                        # label decode y_true
                        if y_true [i] == k:
                            loss[i] -= raw_prediction[i, k]
                        p[k] /= sum_exps  # p_k = y_pred_k = prob of class k
                        # gradient_k = p_k - (y_true == k)
                        gradient[i, k] = p[k] - (y_true[i] == k)

                free(p)
        else:
            with nogil, parallel(num_threads=n_threads):
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    max_value = p[n_classes]  # p[-2]
                    sum_exps = p[n_classes + 1]  # p[-1]
                    loss[i] = log(sum_exps) + max_value

                    for k in range(n_classes):
                        # label decode y_true
                        if y_true [i] == k:
                            loss[i] -= raw_prediction[i, k]
                        p[k] /= sum_exps  # p_k = y_pred_k = prob of class k
                        # gradient_k = (p_k - (y_true == k)) * sw
                        gradient[i, k] = (p[k] - (y_true[i] == k)) * sample_weight[i]

                    loss[i] *= sample_weight[i]

                free(p)

        return np.asarray(loss), np.asarray(gradient)

    def _gradient(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[:, :] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[:, :] gradient,
        int n_threads=1
    ):
        cdef:
            int i, k
            int n_samples = y_true.shape[0]
            int n_classes = raw_prediction.shape[1]
            Y_DTYPE_C sum_exps
            Y_DTYPE_C*  p  # temporary buffer

        if sample_weight is None:
            # inner loop over n_classes
            with nogil, parallel(num_threads=n_threads):
                # Define private buffer variables as each thread might use its
                # own.
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    sum_exps = p[n_classes + 1]  # p[-1]

                    for k in range(n_classes):
                        p[k] /= sum_exps  # p_k = y_pred_k = prob of class k
                        # gradient_k = y_pred_k - (y_true == k)
                        gradient[i, k] = p[k] - (y_true[i] == k)

                free(p)
        else:
            with nogil, parallel(num_threads=n_threads):
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    sum_exps = p[n_classes + 1]  # p[-1]

                    for k in range(n_classes):
                        p[k] /= sum_exps  # p_k = y_pred_k = prob of class k
                        # gradient_k = (p_k - (y_true == k)) * sw
                        gradient[i, k] = (p[k] - (y_true[i] == k)) * sample_weight[i]

                free(p)

        return np.asarray(gradient)

    def _gradient_hessian(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[:, :] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[:, :] gradient,
        G_DTYPE_C[:, :] hessian,
        int n_threads=1
    ):
        cdef:
            int i, k
            int n_samples = y_true.shape[0]
            int n_classes = raw_prediction.shape[1]
            Y_DTYPE_C sum_exps
            Y_DTYPE_C* p  # temporary buffer

        if sample_weight is None:
            # inner loop over n_classes
            with nogil, parallel(num_threads=n_threads):
                # Define private buffer variables as each thread might use its
                # own.
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    sum_exps = p[n_classes + 1]  # p[-1]

                    for k in range(n_classes):
                        p[k] /= sum_exps  # p_k = y_pred_k = prob of class k
                        # hessian_k = p_k * (1 - p_k)
                        # gradient_k = p_k - (y_true == k)
                        gradient[i, k] = p[k] - (y_true[i] == k)
                        hessian[i, k] = p[k] * (1. - p[k])

                free(p)
        else:
            with nogil, parallel(num_threads=n_threads):
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    sum_exps = p[n_classes + 1]  # p[-1]

                    for k in range(n_classes):
                        p[k] /= sum_exps  # p_k = y_pred_k = prob of class k
                        # gradient_k = (p_k - (y_true == k)) * sw
                        # hessian_k = p_k * (1 - p_k) * sw
                        gradient[i, k] = (p[k] - (y_true[i] == k)) * sample_weight[i]
                        hessian[i, k] = (p[k] * (1. - p[k])) * sample_weight[i]

                free(p)

        return np.asarray(gradient), np.asarray(hessian)


    # This method simplifies the implementation of hessp in linear models,
    # i.e. the matrix-vector product of the full hessian, not only of the
    # diagonal (in the classes) approximation as implemented above.
    def _gradient_proba(
        self,
        Y_DTYPE_C[::1] y_true,
        Y_DTYPE_C[:, :] raw_prediction,
        Y_DTYPE_C[::1] sample_weight,
        G_DTYPE_C[:, :] gradient,
        G_DTYPE_C[:, :] proba,
        int n_threads=1
    ):
        cdef:
            int i, k
            int n_samples = y_true.shape[0]
            int n_classes = raw_prediction.shape[1]
            Y_DTYPE_C sum_exps
            Y_DTYPE_C*  p  # temporary buffer

        if sample_weight is None:
            # inner loop over n_classes
            with nogil, parallel(num_threads=n_threads):
                # Define private buffer variables as each thread might use its
                # own.
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    sum_exps = p[n_classes + 1]  # p[-1]

                    for k in range(n_classes):
                        proba[i, k] = p[k] / sum_exps  # y_pred_k = prob of class k
                        # gradient_k = y_pred_k - (y_true == k)
                        gradient[i, k] = proba[i, k] - (y_true[i] == k)

                free(p)
        else:
            with nogil, parallel(num_threads=n_threads):
                p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) * (n_classes + 2))

                for i in prange(n_samples, schedule='static'):
                    sum_exp_minus_max(i, raw_prediction, p)
                    sum_exps = p[n_classes + 1]  # p[-1]

                    for k in range(n_classes):
                        proba[i, k] = p[k] / sum_exps  # y_pred_k = prob of class k
                        # gradient_k = (p_k - (y_true == k)) * sw
                        gradient[i, k] = (proba[i, k] - (y_true[i] == k)) * sample_weight[i]

                free(p)

        return np.asarray(gradient), np.asarray(proba)
