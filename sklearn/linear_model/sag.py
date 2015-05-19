"""Solvers for Ridge and LogisticRegression using SAG algorithm"""

# Authors: Danny Sullivan <dbsullivan23@gmail.com>
#          Tom Dupre la Tour <tom.dupre-la-tour@m4x.org>
#
# Licence: BSD 3 clause

import numpy as np
import scipy.sparse as sp
import warnings

from ..utils import ConvergenceWarning
from ..utils.fixes import astype
from ..utils.seq_dataset import ArrayDataset, CSRDataset
from .sgd_fast import Log, SquaredLoss
from .sag_fast import sag, get_auto_eta

from .stochastic_gradient import SPARSE_INTERCEPT_DECAY


def make_dataset(X, y, sample_weight, random_state):
    # check which type of Sequential Dataset is needed
    y = astype(y, dtype=np.float64, copy=False)
    seed = random_state.randint(0, np.iinfo(np.int32).max)

    if sp.issparse(X):
        dataset = CSRDataset(X.data, X.indptr, X.indices,
                             y, sample_weight, seed=seed)
        intercept_decay = SPARSE_INTERCEPT_DECAY
    else:
        dataset = ArrayDataset(X, y, sample_weight, seed=seed)
        intercept_decay = 1.0

    return dataset, intercept_decay


def sag_ridge(X, y, sample_weight=None, alpha=1e-4, max_iter=1000, tol=0.001,
              verbose=0, random_state=None):
    """SAG solver for Ridge regression

    SAG stands for Stochastic Average Gradient: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a constant learning rate.

    IMPORTANT NOTE: Ridge regressor with 'sag' solver converges faster on
    columns that are on the same scale. You can normalize the data by using
    sklearn.preprocessing.StandardScaler on your data before passing it to the
    fit method.

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using the squared euclidean norm
    L2.

    This implementation works with data represented as dense numpy arrays or
    sparse scipy arrays of floating point values for the features.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Training data

    y : numpy array, shape (n_samples,)
        Target values

    sample_weight : array-like, shape (n_samples,), optional
        Weights applied to individual samples (1. for unweighted).

    alpha : float, optional
        Constant that multiplies the regularization term. Defaults to 0.0001

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
        shuffling the data. Used in 'sag' solver.

    Returns
    -------
    coef_ : array, shape (n_features)
        Weight vector.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = linear_model.Ridge(solver='sag')
    >>> clf.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    Ridge(alpha=1.0, copy_X=True, fit_intercept=True, max_iter=None,
          normalize=False, random_state=None, solver='sag', tol=0.001)

    Reference
    ---------
    Schmidt, M., Roux, N. L., & Bach, F. (2013).
    Minimizing finite sums with the stochastic average gradient
    https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf

    See also
    --------
    Ridge, SGDRegressor, ElasticNet, Lasso, SVR
    """
    if max_iter is None:
        max_iter = 1000

    n_samples, n_features = X.shape[0], X.shape[1]
    alpha = float(alpha) / n_samples
    fit_intercept = False

    # initialization
    if sample_weight is None:
        sample_weight = np.ones(n_samples, dtype=np.float64, order='C')
    weight_pos = 1.0
    weight_neg = 1.0

    # These parameters could be saved between two similar fits, to speed up
    # convergence (see warm starting behavior of sag_logistic()).
    # But current Ridge does not have a warm_start parameter, and current
    # RidgeCV never uses sag solver.
    coef_ = np.zeros(n_features, dtype=np.float64, order='C')
    intercept_init = 0.0
    intercept_sum_gradient_init = 0.0
    sum_gradient_init = np.zeros(n_features, dtype=np.float64, order='C')
    gradient_memory_init = np.zeros(n_samples, dtype=np.float64, order='C')
    seen_init = np.zeros(n_samples, dtype=np.int32, order='C')
    num_seen_init = 0

    dataset, intercept_decay = make_dataset(X, y, sample_weight, random_state)

    # set the step_size at 1 / (alpha + L) where L is the max sum of
    # squares for over all samples
    step_size = get_auto_eta(dataset, alpha, n_samples, SquaredLoss(),
                             fit_intercept)
    if step_size * alpha == 1:
        raise ZeroDivisionError("Current sag implementation does not handle "
                                "the case step_size * alpha == 1")

    intercept_, num_seen, max_iter_reached, intercept_sum_gradient = \
        sag(dataset, coef_.ravel(),
            intercept_init, n_samples,
            n_features, tol,
            max_iter,
            SquaredLoss(),
            step_size, alpha,
            sum_gradient_init.ravel(),
            gradient_memory_init.ravel(),
            seen_init.ravel(),
            num_seen_init,
            weight_pos, weight_neg,
            fit_intercept,
            intercept_sum_gradient_init,
            intercept_decay,
            verbose)

    if max_iter_reached:
        warnings.warn("The max_iter was reached which means "
                      "the coef_ did not converge", ConvergenceWarning)

    return coef_


def sag_logistic(X, y, sample_weight=None, alpha=1e-4, max_iter=1000,
                 tol=0.001, verbose=0, random_state=None,
                 warm_start_mem=dict()):
    """SAG solver for LogisticRegression

    SAG stands for Stochastic Average Gradient: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a constant learning rate.

    IMPORTANT NOTE: LinearRegression with with 'sag' solver converges faster on
    columns that are on the same scale. You can normalize the data by using
    sklearn.preprocessing.StandardScaler on your data before passing it to the
    fit method.

    This implementation works with data represented as dense numpy arrays or
    sparse scipy arrays of floating point values for the features. It will
    fit the data according to log loss.

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using either the squared euclidean norm
    L2.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Training data

    y : numpy array, shape (n_samples,)
        Target values

    sample_weight : array-like, shape (n_samples,), optional
            Weights applied to individual samples (1. for unweighted).

    alpha : float, optional
        Constant that multiplies the regularization term. Defaults to 0.0001

    max_iter: int, optional
        The max number of passes over the training data if the stopping
        criterea is not reached. Defaults to 1000.

    tol: double, optional
        The stopping criterea for the weights. The iterations will stop when
        max(change in weights) / max(weights) < tol. Defaults to 0.001

    verbose: integer, optional
        The verbosity level

    random_state: int or numpy.random.RandomState, optional
        The random_state of the pseudo random number generator to use when
        sampling the data.

    warm_start_mem: dict, optional
        The initialization parameters used for warm starting. It is used for
        example in LogisticRegresionCV. If an intercept needs to be fitted,
        warm_start_mem must contains a 'coef' key, with an array of length
        (n_features + 1).

    Returns
    -------
    warm_start_mem : dict
        Contains a 'coef' key with the fitted result, and eventually the
        fitted intercept at the end of the array. Contains also other keys
        used for warm starting.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> clf = linear_model.LogisticRegression(solver='sag')
    >>> clf.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    LogisticRegression(C=1.0, class_weight=None, dual=False,
        fit_intercept=True, intercept_scaling=1, max_iter=100,
        multi_class='ovr', penalty='l2', random_state=None, solver='sag',
        tol=0.0001, verbose=0, warm_start=False)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    Reference
    ---------
    Schmidt, M., Roux, N. L., & Bach, F. (2013).
    Minimizing finite sums with the stochastic average gradient
    https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf

    See also
    --------
    LogisticRegression, SGDClassifier, LinearSVC, Perceptron

    """
    n_samples, n_features = X.shape[0], X.shape[1]
    alpha = float(alpha) / n_samples

    if sample_weight is None:
        sample_weight = np.ones(n_samples, dtype=np.float64, order='C')

    if 'coef' in warm_start_mem.keys():
        coef_init = warm_start_mem['coef']
    else:
        coef_init = np.zeros(n_features, dtype=np.float64, order='C')

    # coef_init contains eventually the intercept_init at the end.
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

    weight_pos = 1.0
    weight_neg = 1.0

    dataset, intercept_decay = make_dataset(X, y, sample_weight, random_state)

    # set the step_size at 4 / (4 * alpha + L + fit_intercept) where L is the
    # max sum of squares for over all samples
    step_size = get_auto_eta(dataset, alpha, n_samples, Log(),
                             fit_intercept)
    if step_size * alpha == 1.:
        raise ZeroDivisionError("Current sag implementation does not handle "
                                "the case step_size * alpha == 1")

    intercept_, num_seen, max_iter_reached, intercept_sum_gradient = \
        sag(dataset, coef_init.ravel(),
            intercept_init, n_samples,
            n_features, tol,
            max_iter,
            Log(),
            step_size, alpha,
            sum_gradient_init.ravel(),
            gradient_memory_init.ravel(),
            seen_init.ravel(),
            num_seen_init,
            weight_pos, weight_neg,
            fit_intercept,
            intercept_sum_gradient_init,
            intercept_decay,
            verbose)

    if max_iter_reached:
        warnings.warn("The max_iter was reached which means "
                      "the coef_ did not converge", ConvergenceWarning)
    if fit_intercept:
        coef_init = np.append(coef_init, intercept_)

    warm_start_mem = {'coef': coef_init, 'sum_gradient': sum_gradient_init,
                      'intercept_sum_gradient': intercept_sum_gradient,
                      'gradient_memory': gradient_memory_init,
                      'seen': seen_init, 'num_seen': num_seen}

    return warm_start_mem
