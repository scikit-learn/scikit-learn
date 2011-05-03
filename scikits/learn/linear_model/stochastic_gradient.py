# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Implementation of Stochastic Gradient Descent (SGD) with dense data."""

import numpy as np

from ..externals.joblib import Parallel, delayed
from .base import BaseSGDClassifier, BaseSGDRegressor
from .sgd_fast import plain_sgd


class SGDClassifier(BaseSGDClassifier):
    """Linear model fitted by minimizing a regularized empirical loss with SGD.

    SGD stands for Stochastic Gradient Descent: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a decreasing strength schedule (aka learning rate).

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using either the squared euclidean norm
    L2 or the absolute norm L1 or a combination of both (Elastic Net). If the
    parameter update crosses the 0.0 value because of the regularizer, the
    update is truncated to 0.0 to allow for learning sparse models and achieve
    online feature selection.

    This implementation works with data represented as dense numpy arrays of
    floating point values for the features.

    Parameters
    ----------
    loss : str, 'hinge' or 'log' or 'modified_huber'
        The loss function to be used. Defaults to 'hinge'. The hinge loss is a
        margin loss used by standard linear SVM models. The 'log' loss is the
        loss of logistic regression models and can be used for probability
        estimation in binary classifiers. 'modified_huber' is another smooth
        loss that brings tolerance to outliers.

    penalty : str, 'l2' or 'l1' or 'elasticnet'
        The penalty (aka regularization term) to be used. Defaults to 'l2' which
        is the standard regularizer for linear SVM models. 'l1' and 'elasticnet'
        migh bring sparsity to the model (feature selection) not achievable with
        'l2'.

    alpha : float
        Constant that multiplies the regularization term. Defaults to 0.0001

    rho : float
        The Elastic Net mixing parameter, with 0 < rho <= 1.
        Defaults to 0.85.

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter: int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle: bool, optional
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    seed: int, optional
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose: integer, optional
        The verbosity level

    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    Attributes
    ----------
    `coef_` : array, shape = [1, n_features] if n_classes == 2 else [n_classes,
    n_features]
        Weights assigned to the features.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn import linear_model
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> Y = np.array([1, 1, 2, 2])
    >>> clf = linear_model.SGDClassifier()
    >>> clf.fit(X, Y)
    SGDClassifier(loss='hinge', n_jobs=1, shuffle=False, verbose=0, n_iter=5,
           fit_intercept=True, penalty='l2', seed=0, rho=1.0, alpha=0.0001)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    LinearSVC, LogisticRegression

    """

    def _fit_binary(self, X, y):
        """Fit a single binary classifier"""
        # interprete X as dense array
        X = np.asanyarray(X, dtype=np.float64, order='C')

        # encode original class labels as 1 (classes[1]) or -1 (classes[0]).
        y_new = np.ones(y.shape, dtype=np.float64, order='C') * -1.0
        y_new[y == self.classes[1]] = 1.0
        y = y_new

        coef_, intercept_ = plain_sgd(self.coef_,
                                      self.intercept_,
                                      self.loss_function,
                                      self.penalty_type,
                                      self.alpha, self.rho,
                                      X, y,
                                      self.n_iter,
                                      int(self.fit_intercept),
                                      int(self.verbose),
                                      int(self.shuffle),
                                      self.seed,
                                      self.class_weight[1],
                                      self.class_weight[0],
                                      self.sample_weight)

        self.coef_ = np.atleast_2d(coef_)
        self.intercept_ = np.asarray(intercept_)

    def _fit_multiclass(self, X, y):
        """Fit a multi-class classifier by combining binary classifiers

        Each binary classifier predicts one class versus all others. This
        strategy is called OVA: One Versus All.
        """
        X = np.asanyarray(X, dtype=np.float64, order='C')

        # Use joblib to run OVA in parallel.
        res = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                delayed(_train_ova_classifier)(i, c, X, y, self.coef_[i],
                                               self.intercept_[i],
                                               self.loss_function,
                                               self.penalty_type, self.alpha,
                                               self.rho, self.n_iter,
                                               self.fit_intercept,
                                               self.verbose, self.shuffle,
                                               self.seed,
                                               self.class_weight[i],
                                               self.sample_weight)
            for i, c in enumerate(self.classes))

        for i, coef, intercept in res:
            self.coef_[i] = coef
            self.intercept_[i] = intercept

    def decision_function(self, X):
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


def _train_ova_classifier(i, c, X, y, coef_, intercept_, loss_function,
                          penalty_type, alpha, rho, n_iter, fit_intercept,
                          verbose, shuffle, seed, class_weight_pos,
                          sample_weight):
    """Inner loop for One-vs.-All scheme"""
    y_i = np.ones(y.shape, dtype=np.float64, order='C') * -1.0
    y_i[y == c] = 1.0
    coef, intercept = plain_sgd(coef_, intercept_, loss_function,
                                penalty_type, alpha, rho,
                                X, y_i, n_iter, fit_intercept,
                                verbose, shuffle, seed, class_weight_pos, 1.0,
                                sample_weight)
    return (i, coef, intercept)


class SGDRegressor(BaseSGDRegressor):
    """Linear model fitted by minimizing a regularized empirical loss with SGD

    SGD stands for Stochastic Gradient Descent: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a decreasing strength schedule (aka learning rate).

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using either the squared euclidean norm
    L2 or the absolute norm L1 or a combination of both (Elastic Net). If the
    parameter update crosses the 0.0 value because of the regularizer, the
    update is truncated to 0.0 to allow for learning sparse models and achieve
    online feature selection.

    This implementation works with data represented as dense numpy arrays of
    floating point values for the features.

    Parameters
    ----------
    loss : str, 'squared_loss' or 'huber'
        The loss function to be used. Defaults to 'squared_loss' which refers
        to the ordinary least squares fit. 'huber' is an epsilon insensitive loss
        function for robust regression.

    penalty : str, 'l2' or 'l1' or 'elasticnet'
        The penalty (aka regularization term) to be used. Defaults to 'l2' which
        is the standard regularizer for linear SVM models. 'l1' and 'elasticnet'
        migh bring sparsity to the model (feature selection) not achievable with
        'l2'.

    alpha : float
        Constant that multiplies the regularization term. Defaults to 0.0001

    rho : float
        The Elastic Net mixing parameter, with 0 < rho <= 1.
        Defaults to 0.85.

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter: int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle: bool, optional
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    seed: int, optional
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose: integer, optional
        The verbosity level.

    p : float
        Epsilon in the epsilon-insensitive huber loss function;
        only if `loss=='huber'`.

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        Weights asigned to the features.

    `intercept_` : array, shape = [1]
        The intercept term.

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn import linear_model
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = linear_model.SGDRegressor()
    >>> clf.fit(X, y)
    SGDRegressor(loss='squared_loss', shuffle=False, verbose=0, n_iter=5,
           fit_intercept=True, penalty='l2', p=0.1, seed=0, rho=1.0,
           alpha=0.0001)

    See also
    --------
    Ridge, ElasticNet, Lasso, SVR

    """

    def _fit_regressor(self, X, y):
        X = np.asanyarray(X, dtype=np.float64, order='C')
        coef_, intercept_ = plain_sgd(self.coef_,
                                      self.intercept_,
                                      self.loss_function,
                                      self.penalty_type,
                                      self.alpha, self.rho,
                                      X, y,
                                      self.n_iter,
                                      int(self.fit_intercept),
                                      int(self.verbose),
                                      int(self.shuffle),
                                      self.seed,
                                      1.0, 1.0,
                                      self.sample_weight)

        self.coef_ = coef_
        self.intercept_ = np.asarray(intercept_)
