# Authors: Peter Prettenhofer <peter.prettenhofer@gmail.com> (main author)
#          Mathieu Blondel (partial_fit support)
#
# License: BSD Style.
"""Implementation of Stochastic Gradient Descent (SGD) with dense data."""

import numpy as np
import scipy.sparse as sp

from ..externals.joblib import Parallel, delayed

from ..base import RegressorMixin
from ..base import ClassifierMixin
from .base import BaseSGD
from ..utils import atleast2d_or_csr, check_arrays
from ..utils.extmath import safe_sparse_dot
from ..utils import safe_asarray

from .sgd_fast import plain_sgd as plain_sgd_dense
from .sgd_fast_sparse import plain_sgd as plain_sgd_sparse
from .sgd_fast import Hinge, Log, ModifiedHuber, SquaredLoss, Huber


def _tocsr(X):
    """Convert X to CSR matrix, preventing a copy if possible"""
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


class SGDClassifier(BaseSGD, ClassifierMixin):
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
        The loss function to be used. Defaults to 'hinge'. The hinge loss is
        a margin loss used by standard linear SVM models. The 'log' loss is
        the loss of logistic regression models and can be used for
        probability estimation in binary classifiers. 'modified_huber'
        is another smooth loss that brings tolerance to outliers.

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

    learning_rate : string, optional
        The learning rate:
        constant: eta = eta0
        optimal: eta = 1.0/(t+t0) [default]
        invscaling: eta = eta0 / pow(t, power_t)

    eta0 : double
        The initial learning rate [default 0.01].

    power_t : double
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
    >>> Y = np.array([1, 1, 2, 2])
    >>> clf = linear_model.SGDClassifier()
    >>> clf.fit(X, Y)
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
    def __init__(self, loss="hinge", penalty='l2', alpha=0.0001,
                rho=0.85, fit_intercept=True, n_iter=5, shuffle=False,
                verbose=0, n_jobs=1, seed=0, learning_rate="optimal",
                eta0=0.0, power_t=0.5, class_weight=None):
        super(SGDClassifier, self).__init__(loss=loss, penalty=penalty,
                                            alpha=alpha, rho=rho,
                                            fit_intercept=fit_intercept,
                                            n_iter=n_iter, shuffle=shuffle,
                                            verbose=verbose, seed=seed,
                                            learning_rate=learning_rate,
                                            eta0=eta0, power_t=power_t)
        self.class_weight = class_weight
        self.classes = None
        self.n_jobs = int(n_jobs)

    def _set_loss_function(self, loss):
        """Set concrete LossFunction."""
        loss_functions = {
            "hinge": Hinge(),
            "log": Log(),
            "modified_huber": ModifiedHuber(),
        }
        try:
            self.loss_function = loss_functions[loss]
        except KeyError:
            raise ValueError("The loss %s is not supported. " % loss)

    def _set_class_weight(self, class_weight, classes, y):
        """Estimate class weights for unbalanced datasets."""
        if class_weight is None:
            # keep the old class_weight if none provided
            class_weight = self.class_weight
        if class_weight is None or len(class_weight) == 0:
            # uniform class weights
            weight = np.ones(classes.shape[0], dtype=np.float64, order='C')
        elif class_weight == 'auto':
            # proportional to the number of samples in the class
            weight = np.array([1.0 / np.sum(y == i) for i in classes],
                              dtype=np.float64, order='C')
            weight *= classes.shape[0] / np.sum(weight)
        else:
            # user-defined dictionary
            weight = np.ones(classes.shape[0], dtype=np.float64, order='C')
            if not isinstance(class_weight, dict):
                raise ValueError("class_weight must be dict, 'auto', or None,"
                                 " got: %r" % class_weight)
            for c in class_weight:
                i = np.searchsorted(classes, c)
                if classes[i] != c:
                    raise ValueError("Class label %d not present." % c)
                else:
                    weight[i] = class_weight[c]

        self._expanded_class_weight = weight

    def partial_fit(self, X, y, classes=None,
            class_weight=None, sample_weight=None):
        """Fit linear model with Stochastic Gradient Descent.

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]
            Subset of the training data

        y : numpy array of shape [n_samples]
            Subset of the target values

        classes : array, shape = [n_classes]
            Classes across all calls to partial_fit.
            Can be obtained by via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.
            Note that y doesn't need to contain all labels in `classes`.


        class_weight : dict, {class_label : weight} or "auto"
            Weights associated with classes.

            The "auto" mode uses the values of y to automatically adjust
            weights inversely proportional to class frequencies.

            If None, values defined in the previous call to partial_fit
            will be used. If partial_fit was never called before,
            uniform weights are assumed.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples.
            If not provided, uniform weights are assumed.

        Returns
        -------
        self : returns an instance of self.
        """
        X = safe_asarray(X, dtype=np.float64, order="C")
        y = np.asarray(y, dtype=np.float64)

        n_samples, n_features = X.shape
        self._check_fit_data(X, y)

        if self.classes is None and classes is None:
            raise ValueError("classes must be passed on the first call "
                             "to partial_fit.")
        elif classes is not None and self.classes is not None:
            if not np.array_equal(self.classes, np.unique(classes)):
                raise ValueError("`classes` is not the same as on last call "
                                 "to partial_fit.")
        elif classes is not None:
            self.classes = classes

        n_classes = self.classes.shape[0]

        # Allocate datastructures from input arguments
        self._set_class_weight(class_weight, self.classes, y)
        sample_weight = self._validate_sample_weight(sample_weight, n_samples)

        if self.coef_ is None:
            self._allocate_parameter_mem(n_classes, n_features,
                                         coef_init=None, intercept_init=None)

        # delegate to concrete training procedure
        if n_classes > 2:
            self._fit_multiclass(X, y, sample_weight)
        elif n_classes == 2:
            self._fit_binary(X, y, sample_weight)
        else:
            raise ValueError("The number of class labels must be "
                             "greater than one.")

        return self

    def fit(self, X, y, coef_init=None, intercept_init=None,
            class_weight=None, sample_weight=None):
        """Fit linear model with Stochastic Gradient Descent.

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        coef_init : array, shape = [n_classes,n_features]
            The initial coeffients to warm-start the optimization.

        intercept_init : array, shape = [n_classes]
            The initial intercept to warm-start the optimization.

        class_weight : dict, {class_label : weight} or "auto"
            Weights associated with classes. If not given, all classes
            are supposed to have weight one.

            The "auto" mode uses the values of y to automatically adjust
            weights inversely proportional to class frequencies.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples.
            If not provided, uniform weights are assumed.

        Returns
        -------
        self : returns an instance of self.
        """
        X = safe_asarray(X, dtype=np.float64, order="C")
        y = np.asarray(y, dtype=np.float64)

        n_samples, n_features = X.shape
        self._check_fit_data(X, y)

        # sort in asc order; largest class id is positive class
        classes = np.unique(y)
        n_classes = classes.shape[0]

        # Allocate datastructures from input arguments
        self._allocate_parameter_mem(n_classes, n_features,
                                     coef_init, intercept_init)

        self.partial_fit(X, y,
                         classes=classes,
                         sample_weight=sample_weight,
                         class_weight=class_weight)

        # fitting is over, we can now transform coef_ to fortran order
        # for faster predictions
        self._set_coef(self.coef_)

        return self

    def decision_function(self, X):
        """Predict signed 'distance' to the hyperplane (aka confidence score)

        Parameters
        ----------
        X : array, shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] if n_classes == 2 else [n_samples,n_classes]
          The signed 'distances' to the hyperplane(s).
        """
        X = atleast2d_or_csr(X)
        scores = safe_sparse_dot(X, self.coef_.T) + self.intercept_
        if self.classes.shape[0] == 2:
            return np.ravel(scores)
        else:
            return scores

    def predict(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : array or scipy.sparse matrix of shape [n_samples, n_features]
           Whether the numpy.array or scipy.sparse matrix is accepted depends
           on the actual implementation

        Returns
        -------
        array, shape = [n_samples]
           Array containing the predicted class labels.
        """
        scores = self.decision_function(X)
        if self.classes.shape[0] == 2:
            indices = np.array(scores > 0, dtype=np.int)
        else:
            indices = scores.argmax(axis=1)
        return self.classes[np.ravel(indices)]

    def predict_proba(self, X):
        """Predict class membership probability

        Parameters
        ----------
        X : array or scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] if n_classes == 2 else [n_samples,
        n_classes]
            Contains the membership probabilities of the positive class.

        """
        if len(self.classes) != 2:
            raise NotImplementedError("predict_(log_)proba only supported"
                                      " for binary classification")
        elif not isinstance(self.loss_function, Log):
            raise NotImplementedError("predict_(log_)proba only supported when"
                                      " loss='log' (%s given)" % self.loss)

        return 1.0 / (1.0 + np.exp(-self.decision_function(X)))

    def _fit_binary(self, X, y, sample_weight):
        if sp.issparse(X):
            self._fit_binary_sparse(X, y, sample_weight)
        else:
            self._fit_binary_dense(X, y, sample_weight)

    def _fit_multiclass(self, X, y, sample_weight):
        if sp.issparse(X):
            self._fit_multiclass_sparse(X, y, sample_weight)
        else:
            self._fit_multiclass_dense(X, y, sample_weight)

    def _fit_binary_dense(self, X, y, sample_weight):
        """Fit a single binary classifier"""
        # encode original class labels as 1 (classes[1]) or -1 (classes[0]).
        y_new = np.ones(y.shape, dtype=np.float64, order='C') * -1.0
        y_new[y == self.classes[1]] = 1.0
        y = y_new

        coef_, intercept_ = plain_sgd_dense(self.coef_.ravel(),
                                            self.intercept_[0],
                                            self.loss_function,
                                            self.penalty_type,
                                            self.alpha, self.rho,
                                            X, y,
                                            self.n_iter,
                                            int(self.fit_intercept),
                                            int(self.verbose),
                                            int(self.shuffle),
                                            self.seed,
                                            self._expanded_class_weight[1],
                                            self._expanded_class_weight[0],
                                            sample_weight,
                                            self.learning_rate_code,
                                            self.eta0,
                                            self.power_t)

        # need to be 2d
        self.coef_ = coef_.reshape(1, -1)
        # intercept is a float, need to convert it to an array of length 1
        self.intercept_ = np.asarray([intercept_], dtype=np.float64)

    def _fit_multiclass_dense(self, X, y, sample_weight):
        """Fit a multi-class classifier by combining binary classifiers

        Each binary classifier predicts one class versus all others. This
        strategy is called OVA: One Versus All.
        """
        # Use joblib to run OVA in parallel.
        res = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
            delayed(_train_ova_classifier_ds)(i, c, X, y, self.coef_[i],
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

    def _fit_binary_sparse(self, X, y, sample_weight):
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

        coef_, intercept_ = plain_sgd_sparse(self.coef_.ravel(),
                                             self.intercept_[0],
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

        # need to be 2d
        self.coef_ = coef_.reshape(1, -1)
        # intercept is a float, need to convert it to an array of length 1
        self.intercept_ = np.asarray([intercept_], dtype=np.float64)

    def _fit_multiclass_sparse(self, X, y, sample_weight):
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
            delayed(_train_ova_classifier_sp)(i, c, X_data, X_indices,
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


def _train_ova_classifier_ds(i, c, X, y, coef_, intercept_, loss_function,
                             penalty_type, alpha, rho, n_iter, fit_intercept,
                             verbose, shuffle, seed, class_weight_pos,
                             sample_weight, learning_rate, eta0, power_t):
    """Inner loop for one-vs-all scheme."""
    y_i = np.ones(y.shape, dtype=np.float64, order='C') * -1.0
    y_i[y == c] = 1.0
    coef, intercept = plain_sgd_dense(coef_, intercept_, loss_function,
                                      penalty_type, alpha, rho,
                                      X, y_i, n_iter, fit_intercept,
                                      verbose, shuffle, seed,
                                      class_weight_pos, 1.0,
                                      sample_weight, learning_rate, eta0,
                                      power_t)
    return (i, coef, intercept)


def _train_ova_classifier_sp(i, c, X_data, X_indices, X_indptr, y, coef_,
                             intercept_, loss_function, penalty_type, alpha,
                             rho, n_iter, fit_intercept, verbose, shuffle,
                             seed, class_weight_pos, sample_weight,
                             learning_rate, eta0, power_t):
    """Inner loop for One-vs.-All scheme"""
    y_i = np.ones(y.shape, dtype=np.float64, order='C') * -1.0
    y_i[y == c] = 1.0
    coef, intercept = plain_sgd_sparse(coef_, intercept_,
                                       loss_function, penalty_type,
                                       alpha, rho, X_data, X_indices,
                                       X_indptr, y_i, n_iter,
                                       int(fit_intercept), int(verbose),
                                       int(shuffle), int(seed),
                                       class_weight_pos, 1.0,
                                       sample_weight, learning_rate, eta0,
                                       power_t)
    return (i, coef, intercept)


class SGDRegressor(BaseSGD, RegressorMixin):
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
        to the ordinary least squares fit. 'huber' is an epsilon insensitive
        loss function for robust regression.

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
    >>> clf = linear_model.SGDRegressor()
    >>> clf.fit(X, y)
    SGDRegressor(alpha=0.0001, eta0=0.01, fit_intercept=True,
           learning_rate='invscaling', loss='squared_loss', n_iter=5, p=0.1,
           penalty='l2', power_t=0.25, rho=0.85, seed=0, shuffle=False,
           verbose=0)

    See also
    --------
    Ridge, ElasticNet, Lasso, SVR

    """
    def __init__(self, loss="squared_loss", penalty="l2", alpha=0.0001,
                 rho=0.85, fit_intercept=True, n_iter=5, shuffle=False,
                 verbose=0, p=0.1, seed=0, learning_rate="invscaling",
                 eta0=0.01, power_t=0.25):
        self.p = float(p)
        super(SGDRegressor, self).__init__(loss=loss, penalty=penalty,
                                           alpha=alpha, rho=rho,
                                           fit_intercept=fit_intercept,
                                           n_iter=n_iter, shuffle=shuffle,
                                           verbose=verbose, seed=seed,
                                           learning_rate=learning_rate,
                                           eta0=eta0, power_t=power_t)

    def _set_loss_function(self, loss):
        """Get concrete LossFunction"""
        loss_functions = {
            "squared_loss": SquaredLoss(),
            "huber": Huber(self.p),
        }
        try:
            self.loss_function = loss_functions[loss]
        except KeyError:
            raise ValueError("The loss %s is not supported. " % loss)

    def partial_fit(self, X, y, sample_weight=None):
        """Fit linear model with Stochastic Gradient Descent.

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Subset of training data

        y : numpy array of shape [n_samples]
            Subset of target values

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples.
            If not provided, uniform weights are assumed.

        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_arrays(X, y, sparse_format="csr", copy=False)
        y = np.asarray(y, dtype=np.float64, order="C")

        n_samples, n_features = X.shape
        self._check_fit_data(X, y)

        # Allocate datastructures from input arguments
        sample_weight = self._validate_sample_weight(sample_weight, n_samples)

        if self.coef_ is None:
            self._allocate_parameter_mem(1, n_features,
                                         coef_init=None, intercept_init=None)

        self._fit_regressor(X, y, sample_weight)
        return self

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
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

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_arrays(X, y, sparse_format="csr", copy=False)
        y = np.asarray(y, dtype=np.float64, order="C")

        n_samples, n_features = X.shape
        self._check_fit_data(X, y)

        # Allocate datastructures from input arguments
        self._allocate_parameter_mem(1, n_features,
                                     coef_init, intercept_init)

        return self.partial_fit(X, y, sample_weight)

    def decision_function(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : array or scipy.sparse matrix of shape [n_samples, n_features]
           Whether the numpy.array or scipy.sparse matrix is accepted depends
           on the actual implementation.

        Returns
        -------
        array, shape = [n_samples]
           Array containing the predicted class labels.
        """
        X = atleast2d_or_csr(X)
        scores = safe_sparse_dot(X, self.coef_) + self.intercept_
        return scores.ravel()

    def predict(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : array or scipy.sparse matrix of shape [n_samples, n_features]
           Whether the numpy.array or scipy.sparse matrix is accepted depends
           on the actual implementation.

        Returns
        -------
        array, shape = [n_samples]
           Array containing the predicted class labels.
        """
        return self.decision_function(X)

    def _fit_regressor(self, X, y, sample_weight):
        if sp.issparse(X):
            self._fit_regressor_sparse(X, y, sample_weight)
        else:
            self._fit_regressor_dense(X, y, sample_weight)

    def _fit_regressor_dense(self, X, y, sample_weight):
        X = np.asarray(X, dtype=np.float64, order='C')
        coef_, intercept_ = plain_sgd_dense(self.coef_,
                                            self.intercept_[0],
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
                                            sample_weight,
                                            self.learning_rate_code,
                                            self.eta0, self.power_t)

        self.coef_ = coef_
        self.intercept_ = np.asarray([intercept_], dtype=np.float64)

    def _fit_regressor_sparse(self, X, y, sample_weight):
        # interpret X as CSR matrix
        X = _tocsr(X)

        # get sparse matrix datastructures
        X_data = np.array(X.data, dtype=np.float64, order="C")
        X_indices = np.array(X.indices, dtype=np.int32, order="C")
        X_indptr = np.array(X.indptr, dtype=np.int32, order="C")

        coef_, intercept_ = plain_sgd_sparse(self.coef_,
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
        self.intercept_ = np.asarray([intercept_], dtype=np.float64)
