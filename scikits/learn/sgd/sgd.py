# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Implementation of Stochastic Gradient Descent (SGD) with dense data."""

import numpy as np

from ..externals.joblib import Parallel, delayed
from .base import BaseSGD
from .sgd_fast import plain_sgd


class SGD(BaseSGD):
    """Linear model fitted by minimizing a regularized empirical loss with SGD

    SGD stands for Stochastic Gradient Descent: the gradient of the loss is
    estimated each sample at a time and the model is update along the way with a
    decreasing strength schedule (aka learning rate).

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using either the squared euclidean

    This implementation works with data represented as dense numpy arrays of
    floating point values for the features.

    Parameters
    ----------
    loss : str, 'hinge' or 'log' or 'modifiedhuber'
        The loss function to be used. Defaults to 'hinge'.

    penalty : str, 'l2' or 'l1' or 'elasticnet'
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
    >>> from scikits.learn.sgd import SGD
    >>> clf = SGD()
    >>> clf.fit(X, Y)
    SGD(loss='hinge', n_jobs=1, shuffle=False, verbose=0, fit_intercept=True,
      n_iter=5, penalty='l2', coef_=array([ 9.80373,  9.80373]), rho=1.0,
      alpha=0.0001, intercept_=array(-10.0))
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    LinearSVC

    """

    def fit(self, X, y, **params):
        """Fit current model with regularized Stochastic Gradient Descent"""
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float64)
        y = np.asanyarray(y, dtype=np.float64)

        # largest class id is positive class
        classes = np.unique(y)
        self.classes = classes
        if classes.shape[0] > 2:
            self._fit_multiclass(X, y)
        elif classes.shape[0] == 2:
            self._fit_binary(X, y)
        else:
            raise ValueError("The number of class labels must be "
                             "greater than one.")
        # return self for chaining fit and predict calls
        return self

    def _fit_binary(self, X, y):
        """Fit a single binary classifier"""
        # encode original class labels as 1 (classes[1]) or -1 (classes[0]).
        y_new = np.ones(y.shape, dtype=np.float64) * -1.0
        y_new[y == self.classes[1]] = 1.0
        y = y_new

        n_samples, n_features = X.shape[0], X.shape[1]
        if self.coef_ is None:
            self.coef_ = np.zeros(n_features, dtype=np.float64, order="C")
        else:
            if self.coef_.shape != (n_features,):
                raise ValueError("Provided coef_ does not match dataset. ")
        if self.intercept_ is None:
            self.intercept_ = np.zeros(1, dtype=np.float64)
        else:
            if self.intercept_.shape != (1,):
                raise ValueError("Provided intercept_ does not match dataset. ")

        coef_, intercept_ = plain_sgd(self.coef_,
                                      self.intercept_,
                                      self.loss_function,
                                      self.penalty_type,
                                      self.alpha, self.rho,
                                      X, y,
                                      self.n_iter,
                                      int(self.fit_intercept),
                                      int(self.verbose),
                                      int(self.shuffle))

        self.coef_ = coef_
        self.intercept_ = np.asarray(intercept_)

    def _fit_multiclass(self, X, y):
        """Fit a multi-class classifier by combining binary classifiers

        Each binary classifier predicts one class versus all others. This
        strategy is called OVA: One Versus All.
        """
        n_classes = self.classes.shape[0]
        n_samples, n_features = X.shape[0], X.shape[1]
        if self.coef_ is None:
            coef_ = np.zeros((n_classes, n_features),
                             dtype=np.float64, order="C")
        else:
            if self.coef_.shape != (n_classes, n_features):
                raise ValueError("Provided coef_ does not match dataset. ")
            coef_ = self.coef_

        if self.intercept_ is None or isinstance(self.intercept_, float):
            intercept_ = np.zeros(n_classes, dtype=np.float64, order="C")
        else:
            if self.intercept_.shape != (n_classes, ):
                raise ValueError("Provided intercept_ does not match dataset. ")
            intercept_ = self.intercept_

        self.coef_ = coef_
        self.intercept_ = intercept_

        res = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                delayed(_train_ova_classifier)(i, c, X, y, self)
            for i, c in enumerate(self.classes))

        for i, coef, intercept in res:
            coef_[i] = coef
            intercept_[i] = intercept

        self.coef_ = coef_
        self.intercept_ = intercept_

    def predict_margin(self, X):
        """Predict signed 'distance' to the hyperplane (aka confidence score)

        Parameters
        ----------
        X : array, shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] if n_classes == 2 else [n_samples, n_classes]
          The signed 'distances' to the hyperplane(s).
        """
        X = np.atleast_2d(np.asanyarray(X))
        scores = np.dot(X, self.coef_.T) + self.intercept_
        if self.classes.shape[0] == 2:
            return np.ravel(scores)
        else:
            return scores

    def __reduce__(self):
        """Handler which is called at pickeling time

        This is important for joblib because otherwise it will crash trying to
        pickle the external loss function object.
        """
        return SGD, (self.loss, self.penalty, self.alpha, self.rho, self.coef_,
                     self.intercept_, self.fit_intercept, self.n_iter,
                     self.shuffle, self.verbose, self.n_jobs)


def _train_ova_classifier(i, c, X, y, clf):
    """Inner loop for One-vs.-All scheme"""
    y_i = np.ones(y.shape, dtype=np.float64) * -1.0
    y_i[y == c] = 1.0
    coef, intercept = plain_sgd(clf.coef_[i], clf.intercept_[i],
                                clf.loss_function, clf.penalty_type,
                                clf.alpha, clf.rho, X, y_i, clf.n_iter,
                                clf.fit_intercept, clf.verbose,
                                clf.shuffle)
    return (i, coef, intercept)

