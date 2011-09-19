"""
Randomized Lasso: feature selection based on lasso
"""

# Author: Gael Varoquaux
#
# License: BSD Style.

import numpy as np
from scipy.sparse import issparse

from .base import LinearModel
from ..utils import as_float_array, check_random_state, safe_asanyarray
from ..externals.joblib import Parallel, delayed
from .least_angle import lars_path, LassoLarsIC


################################################################################
# Randomized lasso: feature selection
def _randomized_lasso(X, y, weights, alpha, precompute, eps, verbose,
                      max_iter):
    X = (1 + weights) * X
    _, _, coef_ = lars_path(X, y,
                Gram=precompute, overwrite_X=True,
                overwrite_Gram=True, alpha_min=alpha,
                method='lasso', verbose=verbose,
                max_iter=max_iter, eps=eps)
    return coef_[:, -1] != 0

# XXX: must implement this in a more generic way so that it can be
# instantiated for Lasso or LogisticRegression

class RandomizedLasso(LinearModel):
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

    overwrite_X : boolean, optionnal
        If True, X will not be copied
        Default is False

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
        Number of CPUs to use during the cross validation. If '-1', use
        all the CPUs

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

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

    def __init__(self, alpha='aic', a=.2, n_resampling=200,
                 selection_threshold=.5,
                 fit_intercept=True, verbose=False,
                 normalize=True, refit=True, precompute='auto',
                 max_iter=500, overwrite_X=False,
                 eps=np.finfo(np.float).eps, random_state=None,
                 n_jobs=1):
        self.alpha = alpha
        self.a = a
        self.n_resampling = n_resampling
        self.fit_intercept = fit_intercept
        self.max_iter = max_iter
        self.verbose = verbose
        self.normalize = normalize
        self.refit = refit
        self.precompute = precompute
        self.overwrite_X = overwrite_X
        self.eps = eps
        self.random_state = random_state
        self.n_jobs = n_jobs

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

        X = as_float_array(X, overwrite_X=self.overwrite_X)

        X, y, X_mean, y_mean, X_std = self._center_data(X, y,
                                                        self.fit_intercept,
                                                        self.normalize)
        assert self.precompute in (True, False, None, 'auto')
        alpha = self.alpha
        if alpha in ('aic', 'bic'):
            model = LassoLarsIC(precompute=self.precompute,
                                criterion=self.alpha,
                                max_iter=self.max_iter,
                                eps=self.eps)
            model.fit(X, y)
            self.alpha_ = alpha = model.alpha_

        random_state = check_random_state(self.random_state)

        self.scores_ = np.zeros(n_features)
        for active_set in Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                delayed(_randomized_lasso)(X, y,
                        weights=self.a*random_state.random_integers(0,
                                            1, size=(n_features,)),
                        precompute=self.precompute, alpha=alpha,
                        verbose=min(0, self.verbose - 1),
                        max_iter=self.max_iter, eps=self.eps)
                for _ in range(self.n_resampling)):
            self.scores_ += active_set

        self.scores_ /= self.n_resampling

        if self.refit:
            pass
            # XXX: must transform and use OLS to fit
            # XXX: challenge: it is going to interfere with the
            # self.mean_ and self.intercept_
            #self._set_intercept(X_mean, y_mean, X_std)

        return self

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
        return safe_asanyarray(X)[:, self.get_support(indices=issparse(X))]

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
