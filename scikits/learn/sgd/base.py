# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Stochastic Gradient Descent (SGD) with sparse data. """

import numpy as np

from ..base import BaseEstimator

class LinearModel(BaseEstimator):
    """Linear Model trained by minimizing a regularized training
    error using SGD.

    This implementation works on scipy.sparse X and dense coef_.

    Parameters
    ----------
    loss : str, ('hinge'|'log'|'modifiedhuber')
        The loss function to be used. Defaults to 'hinge'. 
    penalty : str, ('l2'|'l1'|'elasticnet')
        The penalty (aka regularization term) to be used. Defaults to 'l2'.
    alpha : float
        Constant that multiplies the regularization term. Defaults to 0.0001
    rho : float
        The Elastic Net mixing parameter, with 0 < rho <= 1.
        Defaults to 0.85.
    coef_ : array, shape = [n_features] if n_classes == 2 else [n_classes, n_features]
        The initial coeffients to warm-start the optimization.
    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        The initial intercept to warm-start the optimization.
    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.
    n_iter: int
        The number of passes over the training data (aka epochs).
        Defaults to 5. 
    shuffle: bool
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.
    verbose: integer, optional
        The verbosity level
    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'.

    Attributes
    ----------
    `coef_` : array, shape = [n_features] if n_classes == 2 else [n_classes, n_features]
        Weights asigned to the features.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> Y = np.array([1, 1, 2, 2])
    >>> from scikits.learn.sgd.sparse import SGD
    >>> clf = SGD()
    >>> clf.fit(X, Y)
    SGD(loss='hinge', n_jobs=1, shuffle=False, verbose=0, fit_intercept=True,
      n_iter=5, penalty='l2', coef_=array([ 9.80373,  9.80373]), rho=1.0,
      alpha=0.0001, intercept_=array(-0.10000000000000001))
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    LinearSVC

    """

    def __init__(self, loss="hinge", penalty='l2', alpha=0.0001,
                 rho=0.85, coef_=None, intercept_=None,
                 fit_intercept=True, n_iter=5, shuffle=False,
                 verbose=0, n_jobs=1):
        self.loss = loss
        self.penalty = penalty
        self.alpha = alpha
        self.rho = rho
        self.coef_ = np.asarray(coef_) if coef_ is not None else None
        self.intercept_ = intercept_
        if self.intercept_ is not None:
            self.intercept_ = np.asarray(intercept_)
        self.fit_intercept = fit_intercept
        self.n_iter = int(n_iter)
        if self.n_iter <= 0:
            raise ValueError("n_iter must be greater than zero.")
        if not isinstance(shuffle, bool):
            raise ValueError("shuffle must be either True or False")
        self.shuffle = shuffle
        self.verbose = verbose
        self.n_jobs = n_jobs
        self._get_loss_function()
        self._get_penalty_type()

    def _get_penalty_type(self):
        penalty_types = {"l2": 2, "l1": 1, "elasticnet": 3}
        try:
            self.penalty_type = penalty_types[self.penalty]
            if self.penalty_type == 2:
                self.rho = 1.0
            elif self.penalty_type == 1:
                self.rho = 0.0
        except KeyError:
            raise ValueError("Penalty %s is not supported. " % self.penalty)

