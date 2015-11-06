"""Solvers for Ridge and LogisticRegression using SAG algorithm"""

# Authors: Tom Dupre la Tour <tom.dupre-la-tour@m4x.org>
#
# Licence: BSD 3 clause

import numpy as np
import warnings

from ..utils import ConvergenceWarning
from ..utils import check_array
from .base import make_dataset
from .sgd_fast import Log, SquaredLoss
from .sag_fast import sag, get_max_squared_sum


def get_auto_step_size(max_squared_sum, alpha_scaled, loss, fit_intercept):
    """Compute automatic step size for SAG solver

    The step size is set to 1 / (alpha_scaled + L + fit_intercept) where L is
    the max sum of squares for over all samples.

    Parameters
    ----------
    max_squared_sum : float
        Maximum squared sum of X over samples.

    alpha_scaled : float
        Constant that multiplies the regularization term, scaled by
        1. / n_samples, the number of samples.

    loss : string, in {"log", "squared"}
        The loss function used in SAG solver.

    fit_intercept : bool
        Specifies if a constant (a.k.a. bias or intercept) will be
        added to the decision function.

    Returns
    -------
    step_size : float
        Step size used in SAG solver.

    """
    if loss == 'log':
        # inverse Lipschitz constant for log loss
        return 4.0 / (max_squared_sum + int(fit_intercept)
                      + 4.0 * alpha_scaled)
    elif loss == 'squared':
        # inverse Lipschitz constant for squared loss
        return 1.0 / (max_squared_sum + int(fit_intercept) + alpha_scaled)
    else:
        raise ValueError("Unknown loss function for SAG solver, got %s "
                         "instead of 'log' or 'squared'" % loss)


def sag_solver(X, y, sample_weight=None, loss='log', alpha=1.,
               max_iter=1000, tol=0.001, verbose=0, random_state=None,
               check_input=True, max_squared_sum=None,
               warm_start_mem=dict()):
    """SAG solver for Ridge and LogisticRegression

    SAG stands for Stochastic Average Gradient: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a constant learning rate.

    IMPORTANT NOTE: 'sag' solver converges faster on columns that are on the
    same scale. You can normalize the data by using
    sklearn.preprocessing.StandardScaler on your data before passing it to the
    fit method.

    This implementation works with data represented as dense numpy arrays or
    sparse scipy arrays of floating point values for the features. It will
    fit the data according to squared loss or log loss.

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using the squared euclidean norm L2.

    .. versionadded:: 0.17

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Training data

    y : numpy array, shape (n_samples,)
        Target values

    sample_weight : array-like, shape (n_samples,), optional
        Weights applied to individual samples (1. for unweighted).

    loss : 'log' | 'squared'
        Loss function that will be optimized.
        'log' is used for classification, like in LogisticRegression.
        'squared' is used for regression, like in Ridge.

    alpha : float, optional
        Constant that multiplies the regularization term. Defaults to 1.

    max_iter: int, optional
        The max number of passes over the training data if the stopping
        criterea is not reached. Defaults to 1000.

    tol: double, optional
        The stopping criterea for the weights. The iterations will stop when
        max(change in weights) / max(weights) < tol. Defaults to .001

    verbose: integer, optional
        The verbosity level.

    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    check_input : bool, default True
        If False, the input arrays X and y will not be checked.

    max_squared_sum : float, default None
        Maximum squared sum of X over samples. If None, it will be computed,
        going through all the samples. The value should be precomputed
        to speed up cross validation.

    warm_start_mem: dict, optional
        The initialization parameters used for warm starting. It is currently
        not used in Ridge.

    Returns
    -------
    coef_ : array, shape (n_features)
        Weight vector.

    n_iter_ : int
        The number of full pass on all samples.

    warm_start_mem : dict
        Contains a 'coef' key with the fitted result, and eventually the
        fitted intercept at the end of the array. Contains also other keys
        used for warm starting.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> X = np.random.randn(n_samples, n_features)
    >>> y = np.random.randn(n_samples)
    >>> clf = linear_model.Ridge(solver='sag')
    >>> clf.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    Ridge(alpha=1.0, copy_X=True, fit_intercept=True, max_iter=None,
          normalize=False, random_state=None, solver='sag', tol=0.001)

    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> clf = linear_model.LogisticRegression(solver='sag')
    >>> clf.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    LogisticRegression(C=1.0, class_weight=None, dual=False,
        fit_intercept=True, intercept_scaling=1, max_iter=100,
        multi_class='ovr', n_jobs=1, penalty='l2', random_state=None,
        solver='sag', tol=0.0001, verbose=0, warm_start=False)

    References
    ----------
    Schmidt, M., Roux, N. L., & Bach, F. (2013).
    Minimizing finite sums with the stochastic average gradient
    https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf

    See also
    --------
    Ridge, SGDRegressor, ElasticNet, Lasso, SVR, and
    LogisticRegression, SGDClassifier, LinearSVC, Perceptron
    """
    # Ridge default max_iter is None
    if max_iter is None:
        max_iter = 1000

    if check_input:
        X = check_array(X, dtype=np.float64, accept_sparse='csr', order='C')
        y = check_array(y, dtype=np.float64, ensure_2d=False, order='C')

    n_samples, n_features = X.shape[0], X.shape[1]
    # As in SGD, the alpha is scaled by n_samples.
    alpha_scaled = float(alpha) / n_samples

    # initialization
    if sample_weight is None:
        sample_weight = np.ones(n_samples, dtype=np.float64, order='C')

    if 'coef' in warm_start_mem.keys():
        coef_init = warm_start_mem['coef']
    else:
        coef_init = np.zeros(n_features, dtype=np.float64, order='C')

    # coef_init contains possibly the intercept_init at the end.
    # Note that Ridge centers the data before fitting, so fit_intercept=False.
    fit_intercept = coef_init.size == (n_features + 1)
    if fit_intercept:
        intercept_init = coef_init[-1]
        coef_init = coef_init[:-1]
    else:
        intercept_init = 0.0

    if 'intercept_sum_gradient' in warm_start_mem.keys():
        intercept_sum_gradient_init = warm_start_mem['intercept_sum_gradient']
    else:
        intercept_sum_gradient_init = 0.0

    if 'gradient_memory' in warm_start_mem.keys():
        gradient_memory_init = warm_start_mem['gradient_memory']
    else:
        gradient_memory_init = np.zeros(n_samples, dtype=np.float64,
                                        order='C')
    if 'sum_gradient' in warm_start_mem.keys():
        sum_gradient_init = warm_start_mem['sum_gradient']
    else:
        sum_gradient_init = np.zeros(n_features, dtype=np.float64, order='C')

    if 'seen' in warm_start_mem.keys():
        seen_init = warm_start_mem['seen']
    else:
        seen_init = np.zeros(n_samples, dtype=np.int32, order='C')

    if 'num_seen' in warm_start_mem.keys():
        num_seen_init = warm_start_mem['num_seen']
    else:
        num_seen_init = 0

    dataset, intercept_decay = make_dataset(X, y, sample_weight, random_state)

    if max_squared_sum is None:
        max_squared_sum = get_max_squared_sum(X)
    step_size = get_auto_step_size(max_squared_sum, alpha_scaled, loss,
                                   fit_intercept)

    if step_size * alpha_scaled == 1:
        raise ZeroDivisionError("Current sag implementation does not handle "
                                "the case step_size * alpha_scaled == 1")

    if loss == 'log':
        class_loss = Log()
    elif loss == 'squared':
        class_loss = SquaredLoss()
    else:
        raise ValueError("Invalid loss parameter: got %r instead of "
                         "one of ('log', 'squared')" % loss)

    intercept_, num_seen, n_iter_, intercept_sum_gradient = \
        sag(dataset, coef_init.ravel(),
            intercept_init, n_samples,
            n_features, tol,
            max_iter,
            class_loss,
            step_size, alpha_scaled,
            sum_gradient_init.ravel(),
            gradient_memory_init.ravel(),
            seen_init.ravel(),
            num_seen_init,
            fit_intercept,
            intercept_sum_gradient_init,
            intercept_decay,
            verbose)

    if n_iter_ == max_iter:
        warnings.warn("The max_iter was reached which means "
                      "the coef_ did not converge", ConvergenceWarning)

    coef_ = coef_init
    if fit_intercept:
        coef_ = np.append(coef_, intercept_)

    warm_start_mem = {'coef': coef_, 'sum_gradient': sum_gradient_init,
                      'intercept_sum_gradient': intercept_sum_gradient,
                      'gradient_memory': gradient_memory_init,
                      'seen': seen_init, 'num_seen': num_seen}

    return coef_, n_iter_, warm_start_mem
