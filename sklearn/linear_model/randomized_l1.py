"""
Randomized Lasso/Logistic: feature selection based on Lasso and
sparse Logistic Regression
"""

# Author: Gael Varoquaux, Alexandre Gramfort
#
# License: BSD 3 clause
import itertools
from abc import ABCMeta, abstractmethod
import warnings

import numpy as np
from scipy.sparse import issparse
from scipy import sparse
from scipy.interpolate import interp1d

from .base import _preprocess_data
from ..base import BaseEstimator, TransformerMixin
from ..externals import six
from ..externals.joblib import Memory, Parallel, delayed
from ..utils import (as_float_array, check_random_state, check_X_y,
                     check_array, safe_mask)
from ..utils.validation import check_is_fitted
from .least_angle import lars_path, LassoLarsIC
from .logistic import LogisticRegression
from ..exceptions import ConvergenceWarning


###############################################################################
# Randomized linear model: feature selection

def _resample_model(estimator_func, X, y, scaling=.5, n_resampling=200,
                    n_jobs=1, verbose=False, pre_dispatch='3*n_jobs',
                    random_state=None, sample_fraction=.75, **params):
    random_state = check_random_state(random_state)
    # We are generating 1 - weights, and not weights
    n_samples, n_features = X.shape

    if not (0 < scaling < 1):
        raise ValueError(
            "'scaling' should be between 0 and 1. Got %r instead." % scaling)

    scaling = 1. - scaling
    scores_ = 0.0
    for active_set in Parallel(n_jobs=n_jobs, verbose=verbose,
                               pre_dispatch=pre_dispatch)(
            delayed(estimator_func)(
                X, y, weights=scaling * random_state.random_integers(
                    0, 1, size=(n_features,)),
                mask=(random_state.rand(n_samples) < sample_fraction),
                verbose=max(0, verbose - 1),
                **params)
            for _ in range(n_resampling)):
        scores_ += active_set

    scores_ /= n_resampling
    return scores_


class BaseRandomizedLinearModel(six.with_metaclass(ABCMeta, BaseEstimator,
                                                   TransformerMixin)):
    """Base class to implement randomized linear models for feature selection

    This implements the strategy by Meinshausen and Buhlman:
    stability selection with randomized sampling, and random re-weighting of
    the penalty.
    """

    @abstractmethod
    def __init__(self):
        pass

    _preprocess_data = staticmethod(_preprocess_data)

    def fit(self, X, y):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training data.

        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        X, y = check_X_y(X, y, ['csr', 'csc'], y_numeric=True,
                         ensure_min_samples=2, estimator=self)
        X = as_float_array(X, copy=False)
        n_samples, n_features = X.shape

        X, y, X_offset, y_offset, X_scale = \
            self._preprocess_data(X, y, self.fit_intercept, self.normalize)

        estimator_func, params = self._make_estimator_and_params(X, y)
        memory = self.memory
        if isinstance(memory, six.string_types):
            memory = Memory(cachedir=memory)

        scores_ = memory.cache(
            _resample_model, ignore=['verbose', 'n_jobs', 'pre_dispatch']
        )(
            estimator_func, X, y,
            scaling=self.scaling, n_resampling=self.n_resampling,
            n_jobs=self.n_jobs, verbose=self.verbose,
            pre_dispatch=self.pre_dispatch, random_state=self.random_state,
            sample_fraction=self.sample_fraction, **params)

        if scores_.ndim == 1:
            scores_ = scores_[:, np.newaxis]
        self.all_scores_ = scores_
        self.scores_ = np.max(self.all_scores_, axis=1)
        return self

    def _make_estimator_and_params(self, X, y):
        """Return the parameters passed to the estimator"""
        raise NotImplementedError

    def get_support(self, indices=False):
        """Return a mask, or list, of the features/indices selected."""
        check_is_fitted(self, 'scores_')

        mask = self.scores_ > self.selection_threshold
        return mask if not indices else np.where(mask)[0]

    # XXX: the two function below are copy/pasted from feature_selection,
    # Should we add an intermediate base class?
    def transform(self, X):
        """Transform a new matrix using the selected features"""
        mask = self.get_support()
        X = check_array(X)
        if len(mask) != X.shape[1]:
            raise ValueError("X has a different shape than during fitting.")
        return check_array(X)[:, safe_mask(X, mask)]

    def inverse_transform(self, X):
        """Transform a new matrix using the selected features"""
        support = self.get_support()
        if X.ndim == 1:
            X = X[None, :]
        Xt = np.zeros((X.shape[0], support.size))
        Xt[:, support] = X
        return Xt


###############################################################################
# Randomized lasso: regression settings

def _randomized_lasso(X, y, weights, mask, alpha=1., verbose=False,
                      precompute=False, eps=np.finfo(np.float).eps,
                      max_iter=500):
    X = X[safe_mask(X, mask)]
    y = y[mask]

    # Center X and y to avoid fit the intercept
    X -= X.mean(axis=0)
    y -= y.mean()

    alpha = np.atleast_1d(np.asarray(alpha, dtype=np.float64))

    X = (1 - weights) * X
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', ConvergenceWarning)
        alphas_, _, coef_ = lars_path(X, y,
                                      Gram=precompute, copy_X=False,
                                      copy_Gram=False, alpha_min=np.min(alpha),
                                      method='lasso', verbose=verbose,
                                      max_iter=max_iter, eps=eps)

    if len(alpha) > 1:
        if len(alphas_) > 1:  # np.min(alpha) < alpha_min
            interpolator = interp1d(alphas_[::-1], coef_[:, ::-1],
                                    bounds_error=False, fill_value=0.)
            scores = (interpolator(alpha) != 0.0)
        else:
            scores = np.zeros((X.shape[1], len(alpha)), dtype=np.bool)
    else:
        scores = coef_[:, -1] != 0.0
    return scores


class RandomizedLasso(BaseRandomizedLinearModel):
    """Randomized Lasso.

    Randomized Lasso works by subsampling the training data and
    computing a Lasso estimate where the penalty of a random subset of
    coefficients has been scaled. By performing this double
    randomization several times, the method assigns high scores to
    features that are repeatedly selected across randomizations. This
    is known as stability selection. In short, features selected more
    often are considered good features.

    Read more in the :ref:`User Guide <randomized_l1>`.

    Parameters
    ----------
    alpha : float, 'aic', or 'bic', optional
        The regularization parameter alpha parameter in the Lasso.
        Warning: this is not the alpha parameter in the stability selection
        article which is scaling.

    scaling : float, optional
        The s parameter used to randomly scale the penalty of different
        features (See :ref:`User Guide <randomized_l1>` for details ).
        Should be between 0 and 1.

    sample_fraction : float, optional
        The fraction of samples to be used in each randomized design.
        Should be between 0 and 1. If 1, all samples are used.

    n_resampling : int, optional
        Number of randomized models.

    selection_threshold: float, optional
        The score above which features should be selected.

    fit_intercept : boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    verbose : boolean or integer, optional
        Sets the verbosity amount

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.
        This parameter is ignored when `fit_intercept` is set to False.
        When the regressors are normalized, note that this makes the
        hyperparameters learned more robust and almost independent of
        the number of samples. The same property is not valid for
        standardized data. However, if you wish to standardize, please
        use `preprocessing.StandardScaler` before calling `fit` on an
        estimator with `normalize=False`.

    precompute : True | False | 'auto'
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter : integer, optional
        Maximum number of iterations to perform in the Lars algorithm.

    eps : float, optional
        The machine-precision regularization in the computation of the
        Cholesky diagonal factors. Increase this for very ill-conditioned
        systems. Unlike the 'tol' parameter in some iterative
        optimization-based algorithms, this parameter does not control
        the tolerance of the optimization.

    n_jobs : integer, optional
        Number of CPUs to use during the resampling. If '-1', use
        all the CPUs

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    memory : Instance of joblib.Memory or string
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    scores_ : array, shape = [n_features]
        Feature scores between 0 and 1.

    all_scores_ : array, shape = [n_features, n_reg_parameter]
        Feature scores between 0 and 1 for all values of the regularization \
        parameter. The reference article suggests ``scores_`` is the max of \
        ``all_scores_``.

    Examples
    --------
    >>> from sklearn.linear_model import RandomizedLasso
    >>> randomized_lasso = RandomizedLasso()

    Notes
    -----
    See examples/linear_model/plot_sparse_recovery.py for an example.

    References
    ----------
    Stability selection
    Nicolai Meinshausen, Peter Buhlmann
    Journal of the Royal Statistical Society: Series B
    Volume 72, Issue 4, pages 417-473, September 2010
    DOI: 10.1111/j.1467-9868.2010.00740.x

    See also
    --------
    RandomizedLogisticRegression, Lasso, ElasticNet
    """
    def __init__(self, alpha='aic', scaling=.5, sample_fraction=.75,
                 n_resampling=200, selection_threshold=.25,
                 fit_intercept=True, verbose=False,
                 normalize=True, precompute='auto',
                 max_iter=500,
                 eps=np.finfo(np.float).eps, random_state=None,
                 n_jobs=1, pre_dispatch='3*n_jobs',
                 memory=Memory(cachedir=None, verbose=0)):
        self.alpha = alpha
        self.scaling = scaling
        self.sample_fraction = sample_fraction
        self.n_resampling = n_resampling
        self.fit_intercept = fit_intercept
        self.max_iter = max_iter
        self.verbose = verbose
        self.normalize = normalize
        self.precompute = precompute
        self.eps = eps
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.selection_threshold = selection_threshold
        self.pre_dispatch = pre_dispatch
        self.memory = memory

    def _make_estimator_and_params(self, X, y):
        assert self.precompute in (True, False, None, 'auto')
        alpha = self.alpha
        if isinstance(alpha, six.string_types) and alpha in ('aic', 'bic'):
            model = LassoLarsIC(precompute=self.precompute,
                                criterion=self.alpha,
                                max_iter=self.max_iter,
                                eps=self.eps)
            model.fit(X, y)
            self.alpha_ = alpha = model.alpha_
        return _randomized_lasso, dict(alpha=alpha, max_iter=self.max_iter,
                                       eps=self.eps,
                                       precompute=self.precompute)


###############################################################################
# Randomized logistic: classification settings

def _randomized_logistic(X, y, weights, mask, C=1., verbose=False,
                         fit_intercept=True, tol=1e-3):
    X = X[safe_mask(X, mask)]
    y = y[mask]
    if issparse(X):
        size = len(weights)
        weight_dia = sparse.dia_matrix((1 - weights, 0), (size, size))
        X = X * weight_dia
    else:
        X *= (1 - weights)

    C = np.atleast_1d(np.asarray(C, dtype=np.float64))
    scores = np.zeros((X.shape[1], len(C)), dtype=np.bool)

    for this_C, this_scores in zip(C, scores.T):
        # XXX : would be great to do it with a warm_start ...
        clf = LogisticRegression(C=this_C, tol=tol, penalty='l1', dual=False,
                                 fit_intercept=fit_intercept)
        clf.fit(X, y)
        this_scores[:] = np.any(
            np.abs(clf.coef_) > 10 * np.finfo(np.float).eps, axis=0)
    return scores


class RandomizedLogisticRegression(BaseRandomizedLinearModel):
    """Randomized Logistic Regression

    Randomized Logistic Regression works by subsampling the training
    data and fitting a L1-penalized LogisticRegression model where the
    penalty of a random subset of coefficients has been scaled. By
    performing this double randomization several times, the method
    assigns high scores to features that are repeatedly selected across
    randomizations. This is known as stability selection. In short,
    features selected more often are considered good features.

    Read more in the :ref:`User Guide <randomized_l1>`.

    Parameters
    ----------
    C : float, optional, default=1
        The regularization parameter C in the LogisticRegression.

    scaling : float, optional, default=0.5
        The s parameter used to randomly scale the penalty of different
        features (See :ref:`User Guide <randomized_l1>` for details ).
        Should be between 0 and 1.

    sample_fraction : float, optional, default=0.75
        The fraction of samples to be used in each randomized design.
        Should be between 0 and 1. If 1, all samples are used.

    n_resampling : int, optional, default=200
        Number of randomized models.

    selection_threshold : float, optional, default=0.25
        The score above which features should be selected.

    fit_intercept : boolean, optional, default=True
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    verbose : boolean or integer, optional
        Sets the verbosity amount

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.
        This parameter is ignored when `fit_intercept` is set to False.
        When the regressors are normalized, note that this makes the
        hyperparameters learnt more robust and almost independent of the number
        of samples. The same property is not valid for standardized data.
        However, if you wish to standardize, please use
        `preprocessing.StandardScaler` before calling `fit` on an estimator
        with `normalize=False`.

    tol : float, optional, default=1e-3
         tolerance for stopping criteria of LogisticRegression

    n_jobs : integer, optional
        Number of CPUs to use during the resampling. If '-1', use
        all the CPUs

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    memory : Instance of joblib.Memory or string
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    scores_ : array, shape = [n_features]
        Feature scores between 0 and 1.

    all_scores_ : array, shape = [n_features, n_reg_parameter]
        Feature scores between 0 and 1 for all values of the regularization \
        parameter. The reference article suggests ``scores_`` is the max \
        of ``all_scores_``.

    Examples
    --------
    >>> from sklearn.linear_model import RandomizedLogisticRegression
    >>> randomized_logistic = RandomizedLogisticRegression()

    Notes
    -----
    See examples/linear_model/plot_sparse_recovery.py for an example.

    References
    ----------
    Stability selection
    Nicolai Meinshausen, Peter Buhlmann
    Journal of the Royal Statistical Society: Series B
    Volume 72, Issue 4, pages 417-473, September 2010
    DOI: 10.1111/j.1467-9868.2010.00740.x

    See also
    --------
    RandomizedLasso, LogisticRegression
    """
    def __init__(self, C=1, scaling=.5, sample_fraction=.75,
                 n_resampling=200,
                 selection_threshold=.25, tol=1e-3,
                 fit_intercept=True, verbose=False,
                 normalize=True,
                 random_state=None,
                 n_jobs=1, pre_dispatch='3*n_jobs',
                 memory=Memory(cachedir=None, verbose=0)):
        self.C = C
        self.scaling = scaling
        self.sample_fraction = sample_fraction
        self.n_resampling = n_resampling
        self.fit_intercept = fit_intercept
        self.verbose = verbose
        self.normalize = normalize
        self.tol = tol
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.selection_threshold = selection_threshold
        self.pre_dispatch = pre_dispatch
        self.memory = memory

    def _make_estimator_and_params(self, X, y):
        params = dict(C=self.C, tol=self.tol,
                      fit_intercept=self.fit_intercept)
        return _randomized_logistic, params

    def _preprocess_data(self, X, y, fit_intercept, normalize=False):
        """Center the data in X but not in y"""
        X, _, X_offset, _, X_scale = _preprocess_data(X, y, fit_intercept,
                                                      normalize=normalize)
        return X, y, X_offset, y, X_scale


###############################################################################
# Stability paths
def _lasso_stability_path(X, y, mask, weights, eps):
    "Inner loop of lasso_stability_path"
    X = X * weights[np.newaxis, :]
    X = X[safe_mask(X, mask), :]
    y = y[mask]

    alpha_max = np.max(np.abs(np.dot(X.T, y))) / X.shape[0]
    alpha_min = eps * alpha_max  # set for early stopping in path
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', ConvergenceWarning)
        alphas, _, coefs = lars_path(X, y, method='lasso', verbose=False,
                                     alpha_min=alpha_min)
    # Scale alpha by alpha_max
    alphas /= alphas[0]
    # Sort alphas in assending order
    alphas = alphas[::-1]
    coefs = coefs[:, ::-1]
    # Get rid of the alphas that are too small
    mask = alphas >= eps
    # We also want to keep the first one: it should be close to the OLS
    # solution
    mask[0] = True
    alphas = alphas[mask]
    coefs = coefs[:, mask]
    return alphas, coefs


def lasso_stability_path(X, y, scaling=0.5, random_state=None,
                         n_resampling=200, n_grid=100,
                         sample_fraction=0.75,
                         eps=4 * np.finfo(np.float).eps, n_jobs=1,
                         verbose=False):
    """Stabiliy path based on randomized Lasso estimates

    Read more in the :ref:`User Guide <randomized_l1>`.

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        training data.

    y : array-like, shape = [n_samples]
        target values.

    scaling : float, optional, default=0.5
        The alpha parameter in the stability selection article used to
        randomly scale the features. Should be between 0 and 1.

    random_state : integer or numpy.random.RandomState, optional
        The generator used to randomize the design.

    n_resampling : int, optional, default=200
        Number of randomized models.

    n_grid : int, optional, default=100
        Number of grid points. The path is linearly reinterpolated
        on a grid between 0 and 1 before computing the scores.

    sample_fraction : float, optional, default=0.75
        The fraction of samples to be used in each randomized design.
        Should be between 0 and 1. If 1, all samples are used.

    eps : float, optional
        Smallest value of alpha / alpha_max considered

    n_jobs : integer, optional
        Number of CPUs to use during the resampling. If '-1', use
        all the CPUs

    verbose : boolean or integer, optional
        Sets the verbosity amount

    Returns
    -------
    alphas_grid : array, shape ~ [n_grid]
        The grid points between 0 and 1: alpha/alpha_max

    scores_path : array, shape = [n_features, n_grid]
        The scores for each feature along the path.

    Notes
    -----
    See examples/linear_model/plot_sparse_recovery.py for an example.
    """
    rng = check_random_state(random_state)

    if not (0 < scaling < 1):
        raise ValueError("Parameter 'scaling' should be between 0 and 1."
                         " Got %r instead." % scaling)

    n_samples, n_features = X.shape

    paths = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_lasso_stability_path)(
            X, y, mask=rng.rand(n_samples) < sample_fraction,
            weights=1. - scaling * rng.random_integers(0, 1,
                                                       size=(n_features,)),
            eps=eps)
        for k in range(n_resampling))

    all_alphas = sorted(list(set(itertools.chain(*[p[0] for p in paths]))))
    # Take approximately n_grid values
    stride = int(max(1, int(len(all_alphas) / float(n_grid))))
    all_alphas = all_alphas[::stride]
    if not all_alphas[-1] == 1:
        all_alphas.append(1.)
    all_alphas = np.array(all_alphas)
    scores_path = np.zeros((n_features, len(all_alphas)))

    for alphas, coefs in paths:
        if alphas[0] != 0:
            alphas = np.r_[0, alphas]
            coefs = np.c_[np.ones((n_features, 1)), coefs]
        if alphas[-1] != all_alphas[-1]:
            alphas = np.r_[alphas, all_alphas[-1]]
            coefs = np.c_[coefs, np.zeros((n_features, 1))]
        scores_path += (interp1d(alphas, coefs,
                        kind='nearest', bounds_error=False,
                        fill_value=0, axis=-1)(all_alphas) != 0)

    scores_path /= n_resampling
    return all_alphas, scores_path
