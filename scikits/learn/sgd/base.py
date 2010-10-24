# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Stochastic Gradient Descent (SGD) with sparse data. """

from ..base import BaseEstimator


class LinearModel(BaseEstimator):
    """Linear Model trained by regularized Stochastic Gradient Descent

    SGD works by iteratively minimizing (sample by sample) the sum a
    running estimate of a loss function (e.g. hinge loss or quadratic
    loss) and a regularizer (e.g. the squared euclidean (L2) norm of the
    coefs) that encodes apriori knowledge of the distribution of the coefs
    (e.g. centered guaussian distribution for the L2 regularizer.

    Parameters
    ----------
    loss : str, ('hinge'|'log'|'modifiedhuber')
        The loss function to be used.
    penalty : str, ('l2'|'l1'|'elasticnet')
        The penalty (aka regularization term) to be used.
    alpha : float
        Constant that multiplies the regularization term. Defaults to 0.0001
    rho : float
        The Elastic Net mixing parameter, with 0 < rho <= 1.
    coef_ : ndarray of shape n_features
        The initial coeffients to warm-start the optimization
    intercept_ : float
        The initial intercept to warm-start the optimization
    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.
    n_iter: int
        The number of passes over the training data (aka epochs).
    shuffle: bool
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        Weights asigned to the features.

    `intercept_` : float
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> Y = np.array([1, 1, 2, 2])
    >>> from scikits.learn.sgd.sparse import SGD
    >>> clf = SGD()
    >>> clf.fit(X, Y)
    SGD(loss='hinge', shuffle=False, fit_intercept=True, n_iter=5, penalty='l2',
      coef_=array([-9.80373, -9.80373]), rho=1.0, alpha=0.0001, intercept_=0.1)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    LinearSVC

    """

    def __init__(self, loss="hinge", penalty='l2', alpha=0.0001,
                 rho=0.85, coef_=None, intercept_=0.0,
                 fit_intercept=True, n_iter=5, shuffle=False):
        self.loss = loss
        self.penalty = penalty
        self.alpha = alpha
        self.rho = rho
        self.coef_ = coef_
        self.intercept_ = intercept_
        self.fit_intercept = fit_intercept
        self.n_iter = int(n_iter)
        if self.n_iter <= 0:
            raise ValueError("n_iter must be greater than zero.")
        if not isinstance(shuffle, bool):
            raise ValueError("shuffle must be either True or False")
        self.shuffle = shuffle
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

