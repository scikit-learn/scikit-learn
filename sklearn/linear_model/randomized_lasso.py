"""
Randomized Lasso: feature selection based on lasso
"""

# Author: Gael Varoquaux
#
# License: BSD Style.

import numpy as np
from scipy.sparse import issparse

from .base import center_data
from ..base import TransformerMixin
from ..utils import as_float_array, check_random_state, safe_asarray
from ..externals.joblib import Parallel, delayed
from .least_angle import lars_path, LassoLarsIC
from .logistic import LogisticRegression
from ..externals.joblib import Memory


################################################################################
# Randomized linear model: feature selection

def _resample_model(estimator_func, X, y, a=.5, n_resampling=200,
                    n_jobs=1, verbose=False, pre_dispatch='3*n_jobs',
                    random_state=None, jacknife_fraction=.75, **params):
    random_state = check_random_state(random_state)
    # We are generating 1 - weights, and not weights
    n_samples, n_features = X.shape
    a = 1 - a
    scores_ = np.zeros(n_features)
    for active_set in Parallel(n_jobs=n_jobs, verbose=verbose,
                               pre_dispatch=pre_dispatch)(
                delayed(estimator_func)(X, y,
                        weights=a * random_state.random_integers(0,
                                                    1, size=(n_features,)),
                        mask=(random_state.rand(n_samples) < jacknife_fraction),
                        verbose=max(0, verbose - 1),
                        **params)
                for _ in range(n_resampling)):
        scores_ += active_set

    scores_ /= n_resampling
    return scores_


class BaseRandomizedLinearModel(TransformerMixin):
    """ Base class to implement randomized linear models for feature
        selection, in the spirit of Meinshausen and Buhlman's:
        stability selection with jacknife, and random reweighting of
        the penalty
    """

    _center_data = staticmethod(center_data)

    def __init__(self, jacknife_fraction=.75, n_resampling=200,
                 selection_threshold=.25, fit_intercept=True, verbose=False,
                 normalize=True, random_state=None, n_jobs=1,
                 pre_dispatch='3*n_jobs',
                 memory=Memory(cachedir=None, verbose=0)):
        self.jacknife_fraction = jacknife_fraction
        self.n_resampling = n_resampling
        self.fit_intercept = fit_intercept
        self.verbose = verbose
        self.normalize = normalize
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.selection_threshold = selection_threshold
        self.pre_dispatch = pre_dispatch
        self.memory = memory

    def fit(self, X, y):
        """Fit the model using X, y as training data.

        parameters
        ----------
        x : array-like, shape = [n_samples, n_features]
            training data.

        y : array-like, shape = [n_samples]
            target values.

        returns
        -------
        self : object
            returns an instance of self.
        """

        X = np.atleast_2d(X)
        n_samples, n_features = X.shape
        y = np.atleast_1d(y)

        X = as_float_array(X, overwrite_X=False)

        X, y, X_mean, y_mean, X_std = self._center_data(X, y,
                                                        self.fit_intercept,
                                                        self.normalize)


        estimator_func, params = self._mk_estimator_and_params(X, y)
        memory = self.memory
        if isinstance(memory, basestring):
            memory = Memory(cachedir=memory)

        self.scores_ = memory.cache(_resample_model)(estimator_func, X, y,
                                    a=self.a,
                                    n_resampling=self.n_resampling,
                                    n_jobs=self.n_jobs,
                                    verbose=self.verbose,
                                    pre_dispatch=self.pre_dispatch,
                                    random_state=self.random_state,
                                    jacknife_fraction=self.jacknife_fraction,
                                    **params)
        return self

    def _mk_estimator_and_params(self, X, y):
        """ Return the parameters passed to the estimator
        """
        raise NotImplementedError

    def get_support(self, indices=False):
        """
        Return a mask, or list, of the features/indices selected.
        """
        mask = self.scores_ > self.selection_threshold
        return mask if not indices else np.where(mask)[0]

    # XXX: the two function below are copy/pasted from feature_selection,
    # Should we add an intermediate base class?
    def transform(self, X):
        """
        Transform a new matrix using the selected features
        """
        return safe_asarray(X)[:, self.get_support(indices=issparse(X))]

    def inverse_transform(self, X):
        """
        Transform a new matrix using the selected features
        """
        support_ = self.get_support()
        if X.ndim == 1:
            X = X[None, :]
        Xt = np.zeros((X.shape[0], support_.size))
        Xt[:, support_] = X
        return Xt



################################################################################
# Randomized lasso: regression settings

def _randomized_lasso(X, y, weights, mask, alpha=1., verbose=False,
                      precompute=False, eps=np.finfo(np.float).eps,
                      max_iter=500):
    # XXX: should we refit the intercept?
    X = X[mask]
    X -= X.mean(axis=0)
    y = y[mask]
    X = (1 - weights) * X
    _, _, coef_ = lars_path(X, y,
                Gram=precompute, overwrite_X=True,
                overwrite_Gram=True, alpha_min=alpha,
                method='lasso', verbose=verbose,
                max_iter=max_iter, eps=eps)
    return coef_[:, -1] != 0


class RandomizedLasso(BaseRandomizedLinearModel):
    """TODO

    Parameters
    ----------
    alpha: float, 'aic', or 'bic'
        The regularization parameter alpha parameter in the Lasso.
        Warning: this is not the alpha parameter in the randomized Lasso
        article

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    verbose : boolean or integer, optional
        Sets the verbosity amount

    normalize : boolean, optional
        If True, the regressors X are normalized

    precompute : True | False | 'auto'
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: integer, optional
        Maximum number of iterations to perform in the Lars algorithm.

    eps: float, optional
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

    pre_dispatch: int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediatly
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    memory : Instance of joblib.Memory or string
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    Attributes
    ----------
    TODO

    Examples
    --------
    TODO

    References
    ----------
    TODO

    See also
    --------
    TODO
    """

    def __init__(self, alpha='aic', a=.5, jacknife_fraction=.75,
                 n_resampling=200,
                 selection_threshold=.25,
                 fit_intercept=True, verbose=False,
                 normalize=True, precompute='auto',
                 max_iter=500,
                 eps=np.finfo(np.float).eps, random_state=None,
                 n_jobs=1, pre_dispatch='3*n_jobs',
                 memory=Memory(cachedir=None, verbose=0)):
        self.alpha = alpha
        self.a = a
        self.jacknife_fraction = jacknife_fraction
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

    def _mk_estimator_and_params(self, X, y):
        assert self.precompute in (True, False, None, 'auto')
        alpha = self.alpha
        if alpha in ('aic', 'bic'):
            model = LassoLarsIC(precompute=self.precompute,
                                criterion=self.alpha,
                                max_iter=self.max_iter,
                                eps=self.eps)
            model.fit(X, y)
            self.alpha_ = alpha = model.alpha_
        return _randomized_lasso, dict(alpha=alpha,
                    max_iter=self.max_iter, eps=self.eps,
                    precompute=self.precompute)


################################################################################
# Randomized logistic: classification settings

def _randomized_logistic(X, y, weights, mask, C=1., verbose=False,
                      fit_intercept=True, tol=1e-3):
    X = X[mask]
    y = y[mask]
    X = (1 - weights) * X
    clf = LogisticRegression(C=C, tol=tol, penalty='l1', dual=False,
                fit_intercept=fit_intercept)
    clf.fit(X, y)
    return np.any(np.abs(clf.coef_) > 10 * np.finfo(np.float).eps, axis=0)


class RandomizedLogistic(BaseRandomizedLinearModel):
    """ TODO
    """

    def __init__(self, C=1, a=.5, jacknife_fraction=.75,
                 n_resampling=200,
                 selection_threshold=.25, tol=1e-3,
                 fit_intercept=True, verbose=False,
                 normalize=True,
                 random_state=None,
                 n_jobs=1, pre_dispatch='3*n_jobs',
                 memory=Memory(cachedir=None, verbose=0)):
        self.C = C
        self.a = a
        self.jacknife_fraction = jacknife_fraction
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

    def _mk_estimator_and_params(self, X, y):
        params = dict(C=self.C, tol=self.tol,
                      fit_intercept=self.fit_intercept)
        return _randomized_logistic, params

    def _center_data(self, X, y, fit_intercept, normalize=False):
        """ Center the data in X but not in y
        """
        X, _, Xmean, _, X_std = center_data(X, y, fit_intercept,
                                                 normalize=normalize)
        return X, y, Xmean, y, X_std
