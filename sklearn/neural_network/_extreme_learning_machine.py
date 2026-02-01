"""Extreme Learning Machine"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABC, abstractmethod
from numbers import Real

import numpy as np
from scipy.special import logsumexp

from sklearn.base import (
    BaseEstimator,
    ClassifierMixin,
    MultiOutputMixin,
    RegressorMixin,
    _fit_context,
)
from sklearn.linear_model import Ridge
from sklearn.neural_network._base import ACTIVATIONS, WEIGHTS
from sklearn.preprocessing import LabelBinarizer
from sklearn.utils import column_or_1d
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.multiclass import check_classification_targets, unique_labels
from sklearn.utils.validation import check_is_fitted, validate_data


class ExtremeLearningBase(BaseEstimator, ABC):
    """Base class for ELM classification and regression.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """

    _parameter_constraints: dict = {
        "hidden_layer_sizes": ["array-like"],
        "activation": [
            StrOptions(
                {
                    "identity",
                    "tanh",
                    "relu",
                    "logistic",
                    "softmax",
                    "softmin",
                    "log_sigmoid",
                    "log_softmax",
                }
            )
        ],
        "weight_init": [
            StrOptions(
                {
                    "zeros",
                    "uniform",
                    "normal",
                    "he_uniform",
                    "lecun_uniform",
                    "glorot_uniform",
                    "he_normal",
                    "lecun_normal",
                    "glorot_normal",
                }
            )
        ],
        "direct_links": ["boolean"],
        "random_state": ["random_state"],
        "ridge_alpha": [Interval(Real, None, None, closed="left"), None],
        "rtol": [Interval(Real, None, None, closed="left"), None],
    }

    @abstractmethod
    def __init__(
        self,
        hidden_layer_sizes=(100,),
        activation="identity",
        weight_init="uniform",
        direct_links=False,
        random_state=None,
        ridge_alpha=None,
        rtol=None,
    ):
        self.hidden_layer_sizes = hidden_layer_sizes
        self.activation = activation
        self.direct_links = direct_links
        self.random_state = random_state
        self.weight_init = weight_init
        self.ridge_alpha = ridge_alpha
        self.rtol = rtol

    def _initialize(self):
        # set all attributes, initialize weights etc. for first call
        # Initialize parameters
        fn = ACTIVATIONS[self.activation]
        self._activation_fn = fn
        hidden_layer_sizes = np.asarray(self.hidden_layer_sizes)

        self._weight_mode = WEIGHTS[self.weight_init]

        # weights shape: (n_layers,)
        # biases shape: (n_layers,)
        self.W_ = []
        self.b_ = []
        rng = self._get_generator(self.random_state)

        self.W_.append(
            self._weight_mode(self.n_features_in_, hidden_layer_sizes[0], rng=rng)
        )
        self.b_.append(self._weight_mode(1, hidden_layer_sizes[0], rng=rng).reshape(-1))
        for i, layer in enumerate(hidden_layer_sizes[1:]):
            # (n_hidden, n_features)
            self.W_.append(
                self._weight_mode(
                    hidden_layer_sizes[i],
                    layer,
                    rng=rng,
                )
            )
            # (n_hidden,)
            self.b_.append(
                self._weight_mode(
                    1,
                    layer,
                    rng=rng,
                ).reshape(-1)
            )

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y):
        self._validate_params()
        # Assumption : X, Y have been pre-processed.
        # X shape: (n_samples, n_features)
        # Y shape: (n_samples, n_classes)
        self._initialize()
        # design matrix shape: (n_samples, sum_hidden+n_features)
        # or (n_samples, sum_hidden)
        D = self._design_matrix(X)

        # beta shape: (sum_hidden+n_features, n_classes)
        # or (sum_hidden, n_classes)

        # If ridge_alpha is None, use direct solve using
        # MoorePenrose Pseudo-Inverse, otherwise use ridge regularized form.
        if self.ridge_alpha is None:
            if self.rtol is None:
                self.coef_ = np.linalg.pinv(D) @ y
            else:
                self.coef_ = np.linalg.pinv(D, rcond=self.rtol) @ y
        else:
            ridge = Ridge(alpha=self.ridge_alpha, fit_intercept=False)
            ridge.fit(D, y)
            self.coef_ = ridge.coef_.T
        if self.coef_.ndim == 2 and self.coef_.shape[1] == 1:
            self.coef_ = self.coef_.ravel()
        return self

    def partial_fit(self, X, y):
        # Moore Penrose Pseudoinverse:
        # D+ = (D.T @ D)^-1 @ D.T
        #
        # Least squares solution:
        # D+ @ y = (D.T @ D)^-1 @ D.T @ y
        #
        # We're persisting the gram (D.T @ D) and moment (D.T @ y) matrices
        # and updating them by adding the gram and moment matrices of
        # each consecutive batch.
        #
        # Consider the summation representation of matrix multiplication:
        # (B.T @ B)_ij = sum_k (B.T_ik @ B_kj)
        #
        # Now suppose the full design matrix is formed by vertically stacking
        # two batches, written as [B_1 | B_2]:
        # [B_1 | B_2].T @ [B_1 | B_2] = B_1.T @ B_1 + B_2.T @ B_2
        #
        # It's evident that it is indeed possible to update normal equations
        # through summation

        self._validate_params()
        # Assumption : X, Y have been pre-processed.
        # X shape: (n_samples, n_features)
        # Y shape: (n_samples, n_classes-1)

        if not hasattr(self, "W_"):
            # initialize params only on first call
            self._initialize()

        # design matrix shape: (n_samples, sum_hidden+n_features)
        # or (n_samples, sum_hidden)
        D = self._design_matrix(X)
        if not hasattr(self, "A"):
            self.A = np.zeros((D.shape[1], D.shape[1]))
            self.B = np.zeros((D.shape[1], y.shape[1]))

        self.A += D.T @ D
        self.B += D.T @ y

        # beta shape: (sum_hidden+n_features, n_classes-1)
        # or (sum_hidden, n_classes-1)

        # If ridge_alpha is None, use direct solve using
        # MoorePenrose Pseudo-Inverse, otherwise use ridge regularized form.

        if self.ridge_alpha is None:
            if self.rtol is None:
                self.coef_ = np.linalg.pinv(self.A) @ self.B
            else:
                self.coef_ = np.linalg.pinv(self.A, rcond=self.rtol) @ self.B
        else:
            reg = np.identity(self.A.shape[0]) * self.ridge_alpha
            self.coef_ = np.linalg.solve(self.A + reg, self.B)
        if self.coef_.ndim == 2 and self.coef_.shape[1] == 1:
            self.coef_ = self.coef_.ravel()
        return self

    def _design_matrix(self, X):
        Hs = []
        H_prev = X
        for W, b in zip(self.W_, self.b_, strict=False):
            Z = H_prev @ W.T + b
            self._activation_fn(Z)
            H_prev = Z
            Hs.append(H_prev)
        if self.direct_links:
            Hs.append(X)
        return np.hstack(Hs)

    def _predict(self, X):
        D = self._design_matrix(X)
        return D @ self.coef_

    def _get_generator(self, random_state):
        return np.random.default_rng(random_state)


class ExtremeLearningClassifier(ClassifierMixin, ExtremeLearningBase):
    """Extreme Learning Machine classifier.

    This model fits a feedforward neural network with fixed random hidden-layer
    parameters and solves for the output weights using linear least squares or
    ridge regression. When direct links are enabled, the model architecture corresponds
    to a Random Vector Functional Link (RVFL) network.

    Parameters
    ----------
    hidden_layer_sizes : array-like of shape (n_layers,), default=(100,)
        The ith element represents the number of neurons in the ith
        hidden layer.

    activation : {'identity', 'tanh', 'relu', 'logistic', 'softmax', 'softmin',\
        'log_sigmoid', 'log_softmax'}, default='identity'
        Activation function for the hidden layer(s).

        - 'identity', no-op activation, useful to implement linear bottleneck,
          returns f(x) = x

        - 'tanh', the hyperbolic tan function,
          returns f(x) = tanh(x).

        - 'relu', the rectified linear unit function,
          returns f(x) = max(0, x)

        - 'logistic', the logistic sigmoid function,
          returns f(x) = 1 / (1 + exp(-x)).

        - 'softmax', the K-way softmax function,
          returns f(x) = exp(x - max(x)) / sum(exp(x - max(x)))

        - 'softmin', softmax of the negated input,
          returns f(x) = exp(-x + min(x)) / sum(exp(-x + min(x)))

        - 'log_sigmoid', the logarithm of the sigmoid function,
          returns f(x) = log(1 / (1 + exp(-x)))

        - 'log_softmax', the logarithm of the standard softmax function,
          returns f(x) = x - log(sum(exp(x)))

    weight_init : {'zeros', 'uniform', 'normal', 'he_uniform',\
        'lecun_uniform', 'glorot_uniform', 'he_normal',\
        'lecun_normal', 'glorot_normal'}, default='uniform'
        Distribution used to initialize the random hidden-layer weights.

        The initialization functions generate weight matrices of shape
        (n_hidden_units, n_features), where values are drawn
        according to the selected scheme.

        - 'zeros', initialize all weights to zero

        - 'uniform', draw weights from a uniform distribution over `[0, 1)`

        - 'normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation 1.

        - 'he_uniform', draw weights from a uniform distribution over
          `[-sqrt(6 / n_features), sqrt(6 / n_features))`

        - 'lecun_uniform', draw weights from a uniform distribution over
          `[-sqrt(3 / n_features), sqrt(3 / n_features))`

        - 'glorot_uniform', draw weights from a uniform distribution over
          `[-sqrt(3 / ((n_features + n_hidden_units) / 2)),
          sqrt(3 / ((n_features + n_hidden_units) / 2))]`

        - 'he_normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation `sqrt(2 / n_features)`.

        - 'lecun_normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation `1 / sqrt(n_features)`.

        - 'glorot_normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation
          `sqrt(1 / ((n_features + n_hidden_units) / 2))`.

    direct_links : bool, default=False
        Whether to connect input layer to output nodes.

        When set to `True`, the original input features are concatenated with the
        hidden-layer activations. This corresponds to the Random Vector Functional Link
        (RVFL) architecture.

    random_state : int, RandomState instance, default=None
        Determines random number generation for weights and bias
        initialization.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    ridge_alpha : float, default=None
        Amount of ridge shrinkage to apply in order to improve
        conditioning during Ridge regression. When set to zero or `None`,
        model uses direct solve using Moore-Penrose Pseudo-Inverse.

    rtol : float, default=None
        Cutoff for small singular values for the Moore-Penrose pseudo-inverse. Only
        applies when `ridge_alpha=None`. When set to `None`, the `numpy.linalg.pinv`
        defaults to `rtol=1e-15`.

    Attributes
    ----------
    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    classes_ : ndarray or list of ndarray of shape (n_classes,)
        Class labels for each output.

    W_ : list of ndarray of shape (n_layers,)
        Weight matrices of the hidden layers. The ith element in the list represents the
        weight matrix corresponding to layer i.

    b_ : list of ndarray of shape (n_layers,)
        Bias vectors of the hidden layers. The ith element in the list represents the
        bias term corresponding to layer i.

    coef_ : ndarray of shape (n_features_out, n_outputs)
        Output weight matrix learned by closed-form solver on the
        design matrix (least-squares solution via Moore-Penrose pseudo-inverse,
        or ridge regression when `ridge_alpha` is set).

    See Also
    --------
    ExtremeLearningRegressor : Extreme Learning Machine regressor.

    Notes
    -----
    Classic Extreme Learning Machine (ELM) and Random Vector Functional Link (RVFL)
    formulations typically use a single hidden layer. This implementation generalizes
    these models to support multiple hidden layers, following recent explorations of
    deep ELM/deep RVFL architectures in the literature.

    The solution for the output weights depends on the conditioning of the design
    matrix. Poorly conditioned design matrices may require regularization
    (`ridge_alpha`) or explicit handling of the pseudo-inverse tolerance (`rtol`).

    The `partial_fit` method accumulates the gram and moment matrices (`D.T @ D` and
    `D.T @ y`) and solves a linear system to update the output weights. This update can
    fail with a singular linear system, especially when each batch contains fewer
    samples than the number of features (e.g., many hidden nodes and/or
    `direct_links=True`), or when the design matrix becomes rank deficient due to
    collinear features. In these cases, regularization can improve numerical stability.

    References
    ----------
    Guang-Bin Huang, Qin-Yu Zhu, and Chee-Kheong Siew. "Extreme learning machine: Theory
    and applications." Neurocomputing, Volume 70, Issues 1-3. 2006. 489-501.

    Migel D. Tissera, and Mark D. McDonnell. "Deep extreme learning machines: supervised
    autoencoding architecture for classification." Neurocomputing, Volume 174, Part A.
    2016. 42-49.

    Examples
    --------
    >>> from sklearn.neural_network import ExtremeLearningClassifier
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = make_classification(n_samples=100, random_state=1)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,
    ...                                                     random_state=1)
    >>> clf = ExtremeLearningClassifier(random_state=1).fit(X_train, y_train)
    >>> clf.predict_proba(X_test[:1])
    array([[0.32046412, 0.67953588]])
    >>> clf.predict(X_test[:5, :])
    array([1, 0, 1, 0, 1])
    >>> clf.score(X_test, y_test)
    0.92
    """

    def __init__(
        self,
        hidden_layer_sizes=(100,),
        activation="identity",
        weight_init="uniform",
        direct_links=False,
        random_state=None,
        ridge_alpha=None,
        rtol=None,
    ):
        super().__init__(
            hidden_layer_sizes=hidden_layer_sizes,
            activation=activation,
            weight_init=weight_init,
            direct_links=direct_links,
            random_state=random_state,
            ridge_alpha=ridge_alpha,
            rtol=rtol,
        )

    def fit(self, X, y):
        """Fit the model to data matrix X and target(s) y.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            The input data.

        y : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            The target values (class labels).

        Returns
        -------
        self : object
            Returns the trained gradient free neural network.
        """
        # shape: (n_samples, n_features)
        X, y = validate_data(self, X, y)
        y = column_or_1d(y, warn=True)
        check_classification_targets(y)

        self.classes_ = unique_labels(y)

        self._label_binarizer = LabelBinarizer()
        # (n_samples, n_classes) or (n_samples, 1) for binary
        y = self._label_binarizer.fit_transform(y)

        if y.shape[1] == 1 and len(self.classes_) == 2:
            y = np.hstack([1 - y, y])

        super().fit(X, y)
        return self

    def partial_fit(self, X, y, classes=None):
        """Update the model with a subset of the given data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input data.

        y : array-like of shape (n_samples,)
            The target values.

        classes : array of shape (n_classes,), default=None
            Classes across all calls to partial_fit.
            Can be obtained via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.

        Returns
        -------
        self : object
            Returns the partially trained gradient free neural network.
        """
        # shape: (n_samples, n_features)
        X, y = validate_data(self, X, y, reset=not hasattr(self, "n_features_in_"))
        y = column_or_1d(y, warn=True)
        check_classification_targets(y)

        if not hasattr(self, "classes_"):
            if classes is None:
                raise TypeError("Classes must not be None for first partial_fit call")
            self.classes_ = classes
            self._label_binarizer = LabelBinarizer()
            self._label_binarizer.fit(self.classes_)

        if not set(unique_labels(y)) <= set(self.classes_):
            raise ValueError(
                "Expected only labels in classes_ = %r, "
                "but got %r." % (list(self.classes_), unique_labels(y))
            )

        # (n_samples, n_classes) or (n_samples, 1) for binary
        y = self._label_binarizer.transform(y)
        if y.shape[1] == 1 and len(self.classes_) == 2:
            y = np.hstack([1 - y, y])

        super().partial_fit(X, y)
        return self

    def predict(self, X):
        """Predict using the extreme learning classifier.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y : ndarray, shape (n_samples,) or (n_samples, n_classes)
            The predicted classes.
        """
        check_is_fitted(self)
        X = validate_data(self, X, reset=False)

        if len(self.classes_) == 1:
            return np.full(X.shape[0], self.classes_[0], dtype=self.classes_.dtype)

        out = self._predict_proba(X, check_input=False)
        return self.classes_[np.argmax(out, axis=1)]

    def predict_proba(self, X):
        """Probability estimates.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y_prob : ndarray of shape (n_samples, n_classes)
            The predicted probability of the sample for each class in the
            model, where classes are ordered as they are in `self.classes_`.
        """
        check_is_fitted(self)
        X = validate_data(self, X, reset=False)
        return self._predict_proba(X, check_input=False)

    def _predict_proba(self, X, check_input=True):
        if check_input:
            X = validate_data(self, X, reset=False)
        out = super()._predict(X)
        return np.exp(out - logsumexp(out, axis=1, keepdims=True))


class ExtremeLearningRegressor(RegressorMixin, MultiOutputMixin, ExtremeLearningBase):
    """Extreme Learning Machine regressor.

    This model fits a feedforward neural network with fixed random hidden-layer
    parameters and solves for the output weights using linear least squares or
    ridge regression. When direct links are enabled, the model architecture corresponds
    to a Random Vector Functional Link (RVFL) network.

    Parameters
    ----------
    hidden_layer_sizes : array-like of shape (n_layers,), default=(100,)
        The ith element represents the number of neurons in the ith
        hidden layer.

    activation : {'identity', 'tanh', 'relu', 'logistic', 'softmax', 'softmin',\
        'log_sigmoid', 'log_softmax'}, default='identity'
        Activation function for the hidden layer(s).

        - 'identity', no-op activation, useful to implement linear bottleneck,
          returns f(x) = x

        - 'tanh', the hyperbolic tan function,
          returns f(x) = tanh(x).

        - 'relu', the rectified linear unit function,
          returns f(x) = max(0, x)

        - 'logistic', the logistic sigmoid function,
          returns f(x) = 1 / (1 + exp(-x)).

        - 'softmax', the K-way softmax function,
          returns f(x) = exp(x - max(x)) / sum(exp(x - max(x)))

        - 'softmin', softmax of the negated input,
          returns f(x) = exp(-x + min(x)) / sum(exp(-x + min(x)))

        - 'log_sigmoid', the logarithm of the sigmoid function,
          returns f(x) = log(1 / (1 + exp(-x)))

        - 'log_softmax', the logarithm of the standard softmax function,
          returns f(x) = x - log(sum(exp(x)))

    weight_init : {'zeros', 'uniform', 'normal', 'he_uniform',\
        'lecun_uniform', 'glorot_uniform', 'he_normal',\
        'lecun_normal', 'glorot_normal'}, default='uniform'
        Distribution used to initialize the random hidden-layer weights.

        The initialization functions generate weight matrices of shape
        (n_hidden_units, n_features), where values are drawn
        according to the selected scheme.

        - 'zeros', initialize all weights to zero

        - 'uniform', draw weights from a uniform distribution over `[0, 1)`

        - 'normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation 1.

        - 'he_uniform', draw weights from a uniform distribution over
          `[-sqrt(6 / n_features), sqrt(6 / n_features))`

        - 'lecun_uniform', draw weights from a uniform distribution over
          `[-sqrt(3 / n_features), sqrt(3 / n_features))`

        - 'glorot_uniform', draw weights from a uniform distribution over
          `[-sqrt(3 / ((n_features + n_hidden_units) / 2)),
          sqrt(3 / ((n_features + n_hidden_units) / 2))]`

        - 'he_normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation `sqrt(2 / n_features)`.

        - 'lecun_normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation `1 / sqrt(n_features)`.

        - 'glorot_normal', draw weights from a normal (Gaussian) distribution
          with mean 0 and standard deviation
          `sqrt(1 / ((n_features + n_hidden_units) / 2))`.

    direct_links : bool, default=False
        Whether to connect input layer to output nodes.

        When set to True, the original input features are concatenated with the
        hidden-layer activations. This corresponds to the Random Vector Functional Link
        (RVFL) architecture.

    random_state : int, RandomState instance, default=None
        Determines random number generation for weights and bias
        initialization.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    ridge_alpha : float, default=None
        Amount of ridge shrinkage to apply in order to improve
        conditioning during Ridge regression. When set to zero or `None`,
        model uses direct solve using Moore-Penrose Pseudo-Inverse.

    rtol : float, default=None
        Cutoff for small singular values for the Moore-Penrose pseudo-inverse. Only
        applies when `ridge_alpha=None`. When set to `None`, the `numpy.linalg.pinv`
        defaults to `rtol=1e-15`.

    Attributes
    ----------
    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    W_ : list of ndarray of shape (n_layers,)
        Weight matrices of the hidden layers. The ith element in the list represents the
        weight matrix corresponding to layer i.

    b_ : list of ndarray of shape (n_layers,)
        Bias vectors of the hidden layers. The ith element in the list represents the
        bias term corresponding to layer i.

    coef_ : ndarray of shape (n_features_out, n_outputs)
        Output weight matrix learned by closed-form solver on the
        design matrix (least-squares solution via Moore-Penrose pseudo-inverse,
        or ridge regression when `ridge_alpha` is set).

    See Also
    --------
    ExtremeLearningClassifier : Extreme Learning Machine classifier.

    Notes
    -----
    Classic Extreme Learning Machine (ELM) and Random Vector Functional Link (RVFL)
    formulations typically use a single hidden layer. This implementation generalizes
    these models to support multiple hidden layers, following recent explorations of
    deep ELM/deep RVFL architectures in the literature.

    The solution for the output weights depends on the conditioning of the design
    matrix. Poorly conditioned design matrices may require regularization
    (`ridge_alpha`) or explicit handling of the pseudo-inverse tolerance (`rtol`).

    The `partial_fit` method accumulates the gram and moment matrices (`D.T @ D` and
    `D.T @ y`) and solves a linear system to update the output weights. This update can
    fail with a singular linear system, especially when each batch contains fewer
    samples than the number of features (e.g., many hidden nodes and/or
    `direct_links=True`), or when the design matrix becomes rank deficient due to
    collinear features. In these cases, regularization can improve numerical stability.

    References
    ----------
    Guang-Bin Huang, Qin-Yu Zhu, and Chee-Kheong Siew. "Extreme learning machine: Theory
    and applications." Neurocomputing, Volume 70, Issues 1-3. 2006. 489-501.

    Migel D. Tissera, and Mark D. McDonnell. "Deep extreme learning machines: supervised
    autoencoding architecture for classification." Neurocomputing, Volume 174, Part A.
    2016. 42-49.

    Examples
    --------
    >>> from sklearn.neural_network import ExtremeLearningRegressor
    >>> from sklearn.datasets import make_regression
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = make_regression(n_samples=200, n_features=20, random_state=1)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y,
    ...                                                     random_state=1)
    >>> regr = ExtremeLearningRegressor(random_state=1)
    >>> regr.fit(X_train, y_train)
    ExtremeLearningRegressor(random_state=1)
    >>> regr.predict(X_test[:2])
    array([  18.36756236, -278.01397422])
    >>> regr.score(X_test, y_test)
    1.0
    """

    def __init__(
        self,
        hidden_layer_sizes=(100,),
        activation="identity",
        weight_init="uniform",
        direct_links=False,
        random_state=None,
        ridge_alpha=None,
        rtol=None,
    ):
        super().__init__(
            hidden_layer_sizes=hidden_layer_sizes,
            activation=activation,
            weight_init=weight_init,
            direct_links=direct_links,
            random_state=random_state,
            ridge_alpha=ridge_alpha,
            rtol=rtol,
        )

    def fit(self, X, y):
        """Fit the model to data matrix X and target(s) y.

        Parameters
        ----------
        X : ndarray or sparse matrix of shape (n_samples, n_features)
            The input data.

        y : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            The target values (real numbers).

        Returns
        -------
        self : object
            Returns the trained gradient free neural network.
        """
        X, y = validate_data(self, X, y, multi_output=True)

        if y.ndim == 2 and y.shape[1] == 1:
            y = column_or_1d(y, warn=True)
        if y.ndim == 1:
            y = y.reshape(-1, 1)
        super().fit(X, y)
        return self

    def partial_fit(self, X, y):
        """Update the model with a subset of the given data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input data.

        y : array-like of shape (n_samples,)
            The target values.

        Returns
        -------
        self : object
            Returns the partially trained gradient free neural network.
        """
        # shape: (n_samples, n_features)
        X, y = validate_data(self, X, y, reset=not hasattr(self, "n_features_in_"))
        y = column_or_1d(y, warn=True)
        if y.ndim == 1:
            y = y.reshape(-1, 1)
        super().partial_fit(X, y)
        return self

    def predict(self, X):
        """Predict using the extreme learning regressor.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y : ndarray of shape (n_samples, n_outputs)
            The predicted values.
        """
        check_is_fitted(self)
        X = validate_data(self, X, reset=False)
        out = super()._predict(X)
        return out
