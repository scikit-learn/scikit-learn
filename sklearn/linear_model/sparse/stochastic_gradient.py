# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Implementation of Stochastic Gradient Descent (SGD) with sparse data."""

import numpy as np
import scipy.sparse as sp

from ...externals.joblib import Parallel, delayed
from ..base import BaseSGDClassifier, BaseSGDRegressor
from ..sgd_fast_sparse import plain_sgd

## TODO add flag for intercept learning rate heuristic
##


def _tocsr(X):
    """Convert X to CSR matrix, preventing a copy if possible"""
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


class SGDClassifier(BaseSGDClassifier):
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

    This implementation works on scipy.sparse X and dense `coef_`.

    Parameters
    ----------
    loss : str, 'hinge' or 'log' or 'modified_huber'
        The loss function to be used. Defaults to 'hinge'. The hinge loss is
        a margin loss used by standard linear SVM models. The 'log' loss is
        the loss of logistic regression models and can be used for probability
        estimation in binary classifiers. 'modified_huber' is another smooth
        loss that brings tolerance to outliers.

    penalty : str, 'l2' or 'l1' or 'elasticnet'
        The penalty (aka regularization term) to be used. Defaults to 'l2'
        which is the standard regularizer for linear SVM models. 'l1' and
        'elasticnet' migh bring sparsity to the model (feature selection)
        not achievable with 'l2'.

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

    seed: int, optional
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose: integer, optional
        The verbosity level

    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    learning_rate : string, optional
        The learning rate:
        constant: eta = eta0
        optimal: eta = 1.0/(t+t0) [default]
        invscaling: eta = eta0 / pow(t, power_t)

    eta0 : double, optional
        The initial learning rate [default 0.01].

    power_t : double, optional
        The exponent for inverse scaling learning rate [default 0.25].

    class_weight : dict, {class_label : weight} or "auto" or None, optional
        Preset for the class_weight fit parameter.

        Weights associated with classes. If not given, all classes
        are supposed to have weight one.

        The "auto" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies.

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
    >>> from sklearn import linear_model
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> clf = linear_model.sparse.SGDClassifier()
    >>> clf.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    SGDClassifier(alpha=0.0001, class_weight=None, eta0=0.0,
            fit_intercept=True, learning_rate='optimal', loss='hinge',
            n_iter=5, n_jobs=1, penalty='l2', power_t=0.5, rho=0.85, seed=0,
            shuffle=False, verbose=0)
    >>> print clf.predict([[-0.8, -1]])
    [1]

    See also
    --------
    LinearSVC, LogisticRegression

    """

    def _fit_binary(self, X, y, sample_weight):
        """Fit a binary classifier."""
        X = _tocsr(X)

        # encode original class labels as 1 (classes[1]) or -1 (classes[0]).
        y_new = np.ones(y.shape, dtype=np.float64, order="C") * -1.0
        y_new[y == self.classes[1]] = 1.0
        y = y_new

        # get sparse matrix datastructures
        X_data = np.asarray(X.data, dtype=np.float64, order="C")
        X_indices = np.asarray(X.indices, dtype=np.int32, order="C")
        X_indptr = np.asarray(X.indptr, dtype=np.int32, order="C")

        coef_, intercept_ = plain_sgd(self.coef_,
                                      self.intercept_,
                                      self.loss_function,
                                      self.penalty_type,
                                      self.alpha, self.rho,
                                      X_data,
                                      X_indices, X_indptr, y,
                                      self.n_iter,
                                      int(self.fit_intercept),
                                      int(self.verbose),
                                      int(self.shuffle),
                                      int(self.seed),
                                      self._expanded_class_weight[1],
                                      self._expanded_class_weight[0],
                                      sample_weight,
                                      self.learning_rate_code,
                                      self.eta0, self.power_t)

        self._set_coef(coef_)
        self.intercept_ = np.asarray(intercept_)

    def _fit_multiclass(self, X, y, sample_weight):
        """Fit a multi-class classifier as a combination of binary classifiers

        Each binary classifier predicts one class versus all others
        (OVA: One Versus All).
        """
        X = _tocsr(X)

        # get sparse matrix datastructures
        X_data = np.asarray(X.data, dtype=np.float64, order="C")
        X_indices = np.asarray(X.indices, dtype=np.int32, order="C")
        X_indptr = np.asarray(X.indptr, dtype=np.int32, order="C")

        res = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                delayed(_train_ova_classifier)(i, c, X_data, X_indices,
                                               X_indptr, y, self.coef_[i],
                                               self.intercept_[i],
                                               self.loss_function,
                                               self.penalty_type, self.alpha,
                                               self.rho, self.n_iter,
                                               self.fit_intercept,
                                               self.verbose, self.shuffle,
                                               self.seed,
                                               self._expanded_class_weight[i],
                                               sample_weight,
                                               self.learning_rate_code,
                                               self.eta0, self.power_t)
            for i, c in enumerate(self.classes))

        for i, coef, intercept in res:
            self.coef_[i] = coef
            self.intercept_[i] = intercept

        self._set_coef(self.coef_)


def _train_ova_classifier(i, c, X_data, X_indices, X_indptr, y, coef_,
                          intercept_, loss_function, penalty_type, alpha,
                          rho, n_iter, fit_intercept, verbose, shuffle,
                          seed, class_weight_pos, sample_weight,
                          learning_rate, eta0, power_t):
    """Inner loop for One-vs.-All scheme"""
    y_i = np.ones(y.shape, dtype=np.float64, order='C') * -1.0
    y_i[y == c] = 1.0
    coef, intercept = plain_sgd(coef_, intercept_,
                                loss_function, penalty_type,
                                alpha, rho, X_data, X_indices,
                                X_indptr, y_i, n_iter,
                                int(fit_intercept), int(verbose),
                                int(shuffle), int(seed),
                                class_weight_pos, 1.0,
                                sample_weight, learning_rate, eta0,
                                power_t)
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
    update is truncated to 0.0 to allow for learning sparse models and
    achieve online feature selection.

    This implementation works with data represented as dense numpy arrays
    of floating point values for the features.

    Parameters
    ----------
    loss : str, 'squared_loss' or 'huber'
        The loss function to be used. Defaults to 'squared_loss' which
        refers to the ordinary least squares fit. 'huber' is an epsilon
        insensitive loss function for robust regression.

    penalty : str, 'l2' or 'l1' or 'elasticnet'
        The penalty (aka regularization term) to be used. Defaults to 'l2'
        which is the standard regularizer for linear SVM models. 'l1' and
        'elasticnet' migh bring sparsity to the model (feature selection)
        not achievable with 'l2'.

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

    seed: int, optional
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose: integer, optional
        The verbosity level

    p : float
        Epsilon in the epsilon insensitive huber loss function;
        only if `loss=='huber'`.

    learning_rate : string, optional
        The learning rate:
        constant: eta = eta0
        optimal: eta = 1.0/(t+t0)
        invscaling: eta = eta0 / pow(t, power_t) [default]

    eta0 : double, optional
        The initial learning rate [default 0.01].

    power_t : double, optional
        The exponent for inverse scaling learning rate [default 0.25].

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        Weights asigned to the features.

    `intercept_` : array, shape = [1]
        The intercept term.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = linear_model.sparse.SGDRegressor()
    >>> clf.fit(X, y)
    SGDRegressor(alpha=0.0001, eta0=0.01, fit_intercept=True,
           learning_rate='invscaling', loss='squared_loss', n_iter=5, p=0.1,
           penalty='l2', power_t=0.25, rho=0.85, seed=0, shuffle=False,
           verbose=0)

    See also
    --------
    RidgeRegression, ElasticNet, Lasso, SVR

    """

    def _fit_regressor(self, X, y, sample_weight):
        # interpret X as CSR matrix
        X = _tocsr(X)

        # get sparse matrix datastructures
        X_data = np.array(X.data, dtype=np.float64, order="C")
        X_indices = np.array(X.indices, dtype=np.int32, order="C")
        X_indptr = np.array(X.indptr, dtype=np.int32, order="C")

        coef_, intercept_ = plain_sgd(self.coef_,
                                      self.intercept_,
                                      self.loss_function,
                                      self.penalty_type,
                                      self.alpha, self.rho,
                                      X_data,
                                      X_indices, X_indptr, y,
                                      self.n_iter,
                                      int(self.fit_intercept),
                                      int(self.verbose),
                                      int(self.shuffle),
                                      int(self.seed),
                                      1.0, 1.0,
                                      sample_weight,
                                      self.learning_rate_code,
                                      self.eta0, self.power_t)

        self.coef_ = coef_
        self.intercept_ = np.asarray(intercept_)
