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
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    Attributes
    ----------
    `coef_` : array, shape = [n_features] if n_classes == 2 else [n_classes,
    n_features]
        Weights asigned to the features.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> Y = np.array([1, 1, 2, 2])
    >>> clf = SGDClassifier()
    >>> clf.fit(X, Y)
    SGDClassifier(loss='hinge', n_jobs=1, shuffle=False, verbose=0, n_iter=5,
           fit_intercept=True, penalty='l2', rho=1.0, alpha=0.0001)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    LinearSVC, Logistic

    """

    def fit(self, X, y, coef_init=None, intercept_init=None,
            class_weight={}, **params):
        """Fit linear model with Stochastic Gradient Descent.

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data
        y : numpy array of shape [n_samples]
            Target values
        coef_init : array, shape = [n_features] if n_classes == 2 else [n_classes,
        n_features]
            The initial coeffients to warm-start the optimization.
        intercept_init : array, shape = [1] if n_classes == 2 else [n_classes]
            The initial intercept to warm-start the optimization.
        class_weight : dict, {class_label : weight} or "auto"
            Weights associated with classes. If not given, all classes
            are supposed to have weight one.

            The "auto" mode uses the values of y to automatically adjust
            weights inversely proportional to class frequencies.

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float64)
        y = np.asanyarray(y, dtype=np.float64)

        # largest class id is positive class
        self.classes = np.unique(y)

        self.weight = self._get_class_weight(class_weight, self.classes, y)

        if self.classes.shape[0] > 2:
            self._fit_multiclass(X, y, coef_init, intercept_init)
        elif self.classes.shape[0] == 2:
            self._fit_binary(X, y, coef_init, intercept_init)
        else:
            raise ValueError("The number of class labels must be "
                             "greater than one.")
        # return self for chaining fit and predict calls
        return self

    def _fit_binary(self, X, y, coef_init, intercept_init):
        """Fit a single binary classifier"""
        # encode original class labels as 1 (classes[1]) or -1 (classes[0]).
        y_new = np.ones(y.shape, dtype=np.float64) * -1.0
        y_new[y == self.classes[1]] = 1.0
        y = y_new

        n_samples, n_features = X.shape[0], X.shape[1]
        self.coef_ = np.zeros(n_features, dtype=np.float64, order="C")
        self.intercept_ = np.zeros(1, dtype=np.float64, order="C")
        if coef_init is not None:
            coef_init = np.asanyarray(coef_init)
            if coef_init.shape != (n_features,):
                raise ValueError("Provided coef_init does not match dataset.")
            self.coef_ = coef_init
        if intercept_init is not None:
            intercept_init = np.asanyarray(intercept_init)
            if intercept_init.shape != (1,):
                raise ValueError("Provided intercept_init " \
                                 "does not match dataset.")
            else:
                self.intercept_ = intercept_init

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
                                      self.weight[1], self.weight[0])

        self.coef_ = coef_
        self.intercept_ = np.asarray(intercept_)

    def _fit_multiclass(self, X, y, coef_init, intercept_init):
        """Fit a multi-class classifier by combining binary classifiers

        Each binary classifier predicts one class versus all others. This
        strategy is called OVA: One Versus All.
        """
        n_classes = self.classes.shape[0]
        n_samples, n_features = X.shape[0], X.shape[1]
        self.coef_ = np.zeros((n_classes, n_features),
                         dtype=np.float64, order="C")
        self.intercept_ = np.zeros(n_classes, dtype=np.float64, order="C")

        if coef_init is not None:
            coef_init = np.asanyarray(coef_init)
            if coef_init.shape != (n_classes, n_features):
                raise ValueError("Provided coef_ does not match dataset. ")
            else:
                self.coef_ = coef_init
        if intercept_init is not None:
            intercept_init = np.asanyarray(intercept_init)
            if intercept_init.shape != (n_classes, ):
                raise ValueError("Provided intercept_init " \
                                 "does not match dataset.")
            else:
                self.intercept_ = intercept_init

        res = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                delayed(_train_ova_classifier)(i, c, X, y, self.coef_[i],
                                               self.intercept_[i],
                                               self.loss_function,
                                               self.penalty_type, self.alpha,
                                               self.rho, self.n_iter,
                                               self.fit_intercept,
                                               self.verbose, self.shuffle,
                                               self.weight[i])
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
                          verbose, shuffle, weight_pos):
    """Inner loop for One-vs.-All scheme"""
    y_i = np.ones(y.shape, dtype=np.float64) * -1.0
    y_i[y == c] = 1.0
    coef, intercept = plain_sgd(coef_, intercept_, loss_function,
                                penalty_type, alpha, rho,
                                X, y_i, n_iter, fit_intercept,
                                verbose, shuffle, weight_pos, 1.0)
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

    n_iter: int
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle: bool
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    verbose: integer, optional
        The verbosity level

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        Weights asigned to the features.

    `intercept_` : array, shape = [1]
        The intercept term.

    Examples
    --------
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = SGDRegressor()
    >>> clf.fit(X, y)
    SGDRegressor(loss='squared_loss', shuffle=False, verbose=0, n_iter=5,
           epsilon=0.1, fit_intercept=True, penalty='l2', rho=1.0,
           alpha=0.0001)

    See also
    --------
    RidgeRegression, ElasticNet, Lasso, SVR

    """

    def fit(self, X, y, coef_init=None, intercept_init=None,
            **params):
        """Fit linear model with Stochastic Gradient Descent.

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data
        y : numpy array of shape [n_samples]
            Target values
        coef_init : array, shape = [n_features]
            The initial coeffients to warm-start the optimization.
        intercept_init : array, shape = [1]
            The initial intercept to warm-start the optimization.

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float64)
        y = np.asanyarray(y, dtype=np.float64)

        n_samples, n_features = X.shape[0], X.shape[1]
        self.coef_ = np.zeros(n_features, dtype=np.float64, order="C")
        self.intercept_ = np.zeros(1, dtype=np.float64, order="C")
        if coef_init is not None:
            coef_init = np.asanyarray(coef_init)
            if coef_init.shape != (n_features,):
                raise ValueError("Provided coef_init does not match dataset.")
            self.coef_ = coef_init
        if intercept_init is not None:
            intercept_init = np.asanyarray(intercept_init)
            if intercept_init.shape != (1,):
                raise ValueError("Provided intercept_init " \
                                 "does not match dataset.")
            else:
                self.intercept_ = intercept_init

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
                                      1.0, 1.0)

        self.coef_ = coef_
        self.intercept_ = np.asarray(intercept_)
        return self
