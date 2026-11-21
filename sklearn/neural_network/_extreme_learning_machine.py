"""Extreme Learning Machine"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABC, abstractmethod
from numbers import Real

import numpy as np

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
from sklearn.utils._array_api import (
    check_same_namespace,
    get_namespace,
    get_namespace_and_device,
    move_to,
)
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

    def _initialize(self, X):
        # set all attributes, initialize weights etc. for first call
        # Initialize parameters
        xp, _, device_ = get_namespace_and_device(X)
        fn = ACTIVATIONS[self.activation]
        self._activation_fn = fn
        hidden_layer_sizes = np.asarray(self.hidden_layer_sizes)

        self._weight_mode = WEIGHTS[self.weight_init]

        # weights shape: (n_layers,)
        # biases shape: (n_layers,)
        self.W_ = []
        self.b_ = []
        rng = self._get_generator(self.random_state)

        in_sizes = np.concatenate(([self.n_features_in_], hidden_layer_sizes[:-1]))
        for n_in, n_out in zip(in_sizes, hidden_layer_sizes):
            # (n_hidden, n_features)
            self.W_.append(
                xp.asarray(
                    self._weight_mode(n_in, n_out, rng=rng),
                    dtype=X.dtype,
                    device=device_,
                )
            )
            # (n_hidden,)
            self.b_.append(
                xp.asarray(
                    self._weight_mode(1, n_out, rng=rng).reshape(-1),
                    dtype=X.dtype,
                    device=device_,
                )
            )

    def fit(self, X, y):
        xp, _, _ = get_namespace_and_device(X)

        # Assumption : X, Y have been pre-processed.
        # X shape: (n_samples, n_features)
        # Y shape: (n_samples, n_classes)
        self._initialize(X)
        # design matrix shape: (n_samples, sum_hidden+n_features)
        # or (n_samples, sum_hidden)
        D = self._design_matrix(X)

        # beta shape: (sum_hidden+n_features, n_classes)
        # or (sum_hidden, n_classes)

        # If ridge_alpha is None, use direct solve using
        # MoorePenrose Pseudo-Inverse, otherwise use ridge regularized form.
        if self.ridge_alpha is None:
            self.coef_ = xp.linalg.pinv(D, rtol=self.rtol) @ y
        else:
            ridge = Ridge(alpha=self.ridge_alpha, solver="svd", fit_intercept=False)
            ridge.fit(D, y)
            if ridge.coef_.ndim == 1:
                self.coef_ = ridge.coef_
            else:
                self.coef_ = ridge.coef_.T
        if self.coef_.ndim == 2 and self.coef_.shape[1] == 1:
            self.coef_ = xp.squeeze(self.coef_, axis=1)
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

        xp, _, device_ = get_namespace_and_device(X)

        # Assumption : X, Y have been pre-processed.
        # X shape: (n_samples, n_features)
        # Y shape: (n_samples, n_classes-1)

        if not hasattr(self, "W_"):
            # initialize params only on first call
            self._initialize(X)

        # design matrix shape: (n_samples, sum_hidden+n_features)
        # or (n_samples, sum_hidden)
        D = self._design_matrix(X)
        if not hasattr(self, "A"):
            self.A = xp.zeros((D.shape[1], D.shape[1]), dtype=X.dtype, device=device_)
            self.B = xp.zeros((D.shape[1], y.shape[1]), dtype=X.dtype, device=device_)

        self.A += D.T @ D
        self.B += D.T @ y

        # beta shape: (sum_hidden+n_features, n_classes-1)
        # or (sum_hidden, n_classes-1)

        # If ridge_alpha is None, use direct solve using
        # MoorePenrose Pseudo-Inverse, otherwise use ridge regularized form.

        # A = D.T @ D
        # B = D.T @ y
        # pinv(D) @ y = pinv(A)@B:
        # pinv(D) @ y = (D.T @ D)^-1 @ D.T @ y

        # pinv(A) @ B = pinv(D.T @ D) @ D.T @ y

        # pinv(D.T @ D) @ D.T @ y = ((D.T @ D).T @ (D.T @ D))^-1 @ (D.T @ D).T @ D.T @ y
        # = ((D.T @ D) @ (D.T @ D))^-1 @ (D.T @ D) @ D.T @ y
        # = (D.T @ D)^-1 @ (D.T @ D)^-1 @ (D.T @ D) @ D.T @ y
        # = (D.T @ D)^-1 @ D.T @ y

        if self.ridge_alpha is None:
            self.coef_ = xp.linalg.pinv(self.A, rtol=self.rtol) @ self.B
        else:
            # xp.linalg.solve(self.A + reg_mat, self.B)
            # is equivalent to
            # ridge = Ridge(
            #   alpha=self.ridge_alpha, fit_intercept=False, solver='cholesky'
            #   )
            # ridge.fit(D_all[:current_batch], Y_all[:current_batch])
            ridge_mat = (
                xp.eye(self.A.shape[0], dtype=self.A.dtype, device=device_)
                * self.ridge_alpha
            )

            self.coef_ = xp.linalg.solve(self.A + ridge_mat, self.B)

        self.coef_ = xp.astype(self.coef_, X.dtype)
        if self.coef_.ndim == 2 and self.coef_.shape[1] == 1:
            self.coef_ = xp.squeeze(self.coef_, axis=1)
        return self

    def _design_matrix(self, X):
        xp, _ = get_namespace(X)
        Hs = []
        H_prev = X
        for W, b in zip(self.W_, self.b_, strict=False):
            Z = H_prev @ W.T + b
            self._activation_fn(Z)
            H_prev = Z
            Hs.append(H_prev)
        if self.direct_links:
            Hs.append(X)
        return xp.concat(Hs, axis=1)

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

    @_fit_context(prefer_skip_nested_validation=True)
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
        xp, _, device_ = get_namespace_and_device(X)

        X, y = validate_data(
            self,
            X,
            y,
            dtype=[xp.float64, xp.float32],
        )

        y = column_or_1d(y, warn=True)
        check_classification_targets(y)

        self._label_binarizer = LabelBinarizer(sparse_output=False)
        y = self._label_binarizer.fit_transform(y)
        y = move_to(y, xp=xp, device=device_)
        y = xp.astype(y, X.dtype, copy=False)
        self.classes_ = self._label_binarizer.classes_
        if self.classes_.shape[0] < 2:
            raise ValueError("y must contain more than 1 class.")

        # Preserve two output columns for binary classification.
        if self.classes_.shape[0] == 2:
            y = xp.concat((1 - y, y), axis=1)

        super().fit(X, y)
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def partial_fit(self, X, y, classes=None):
        """Update the model with a subset of the given data.

        Parameters
        ----------

        X : array-like of shape (n_samples, n_features)
          The batched training input samples.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
          The batched target values (class labels).

        classes : array-like of shape (n_classes,), default=None
            Classes across all calls to partial_fit.
            Can be obtained via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.

        Returns
        -------
        self : object
            The partially fitted estimator.

        Notes
        -----
        The design matrix is incrementally updated by persisting the gram matrix
        (D.T @ D) and the moment vector (D.T @ y) for each batch, then adding
        the gram and moment contributions of each new batch.

        For batches D = [D1, D2, ..., Dk].T::

            [D1; D2; ...; Dk].T @ [D1; D2; ...; Dk] = sum_k(Dk.T @ Dk)
            [D1; D2; ...; Dk].T @ [y1; y2; ...; yk] = sum_k(Dk.T @ yk)

        This allows incremental accumulation.

        The way we accumulate information across batches prevents us from using
        ``Ridge()`` as we do in full fit. In ``partial_fit()`` we never have access to
        the full design matrix, we only persist the aggregate quantities
        ``D.T @ D``, ``D.T @ y``. Instead, for the regularized case, we directly solve
        the system using ``scipy.linalg.solve(A + reg_mat, B)``, which is mathematically
        equivalent to scikit-learn's ``Ridge(solver='cholesky')``.

        One other difference between full fit and partial fit arises in the direct
        solve path. Because ``pinv()`` is acting on ``D.T @ D`` as opposed to just D,
        the condition number is squared. This may require a lower ``rtol`` to maintain
        numerical stability and avoid loss of information.
        """
        # shape: (n_samples, n_features)
        xp, _, device_ = get_namespace_and_device(X)

        first_call = not hasattr(self, "classes_")
        if first_call:
            if classes is None:
                raise TypeError("Classes must not be None for first partial_fit call")
            dtype = [xp.float64, xp.float32]

        else:
            xp_W, _, device_W = get_namespace_and_device(self.W_[0])
            if (xp != xp_W) or (device_ != device_W):
                raise ValueError(
                    f"All batches must use same device and namespace , "
                    f"{(xp_W, device_W)!r}, but got {(xp, device_)!r}."
                )
            dtype = self.W_[0].dtype

        X, y = validate_data(
            self,
            X,
            y,
            dtype=dtype,
            reset=first_call,
        )

        y = column_or_1d(y, warn=True)
        check_classification_targets(y)

        if first_call:
            self._label_binarizer = LabelBinarizer(sparse_output=False).fit(classes)
            self.classes_ = self._label_binarizer.classes_
            if self.classes_.shape[0] < 2:
                raise ValueError("classes must contain more than 1 class")

        # Check membership without Python sets or converting labels to NumPy.
        xp_classes, _, device_classes = get_namespace_and_device(self.classes_)
        y = move_to(y, xp=xp_classes, device=device_classes)
        matches = xp_classes.expand_dims(y, axis=1) == xp_classes.expand_dims(
            self.classes_, axis=0
        )
        if not xp_classes.all(xp_classes.any(matches, axis=1)):
            raise ValueError(
                f"Expected only labels in classes_ = {self.classes_!r}, "
                f"but got {unique_labels(y)!r}."
            )

        y = self._label_binarizer.transform(y)
        y = move_to(y, xp=xp, device=device_)
        y = xp.astype(y, X.dtype, copy=False)

        # Preserve two output columns for binary classification.
        if self.classes_.shape[0] == 2:
            y = xp.concat((1 - y, y), axis=1)
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
        check_same_namespace(X, self, attribute="coef_", method="predict")
        xp, _, device_ = get_namespace_and_device(X)

        X = validate_data(
            self,
            X,
            reset=False,
            dtype=[xp.float64, xp.float32],
        )

        out = self._predict_proba(X)
        indices = xp.argmax(out, axis=1)
        xp_classes, _, device_classes = get_namespace_and_device(self.classes_)
        indices = move_to(indices, xp=xp_classes, device=device_classes)
        y_pred = xp_classes.take(self.classes_, indices, axis=0)

        if isinstance(y_pred[0], str):
            return y_pred
        return move_to(y_pred, xp=xp, device=device_)

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
        check_same_namespace(X, self, attribute="coef_", method="predict_proba")
        xp, _ = get_namespace(X)

        X = validate_data(
            self,
            X,
            reset=False,
            dtype=[xp.float64, xp.float32],
        )
        return self._predict_proba(X)

    def _predict_proba(self, X):
        xp, _ = get_namespace(X)
        out = super()._predict(X)
        shifted = out - xp.max(out, axis=1, keepdims=True)
        exp_out = xp.exp(shifted)
        return exp_out / xp.sum(exp_out, axis=1, keepdims=True)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.array_api_support = True
        return tags


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

    @_fit_context(prefer_skip_nested_validation=True)
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
        xp, _, device_ = get_namespace_and_device(X)

        X, y = validate_data(
            self,
            X,
            y,
            dtype=[xp.float64, xp.float32],
            multi_output=True,
            y_numeric=True,
        )

        y = move_to(y, xp=xp, device=device_)
        y = xp.astype(y, X.dtype, copy=False)

        if y.ndim == 1:
            y = xp.expand_dims(y, axis=-1)
        super().fit(X, y)
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def partial_fit(self, X, y):
        """Update the model with a subset of the given data.

        Parameters
        ----------

        X : array-like of shape (n_samples, n_features)
          The batched training input samples.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
          The batched target values (class labels).

        Returns
        -------
        self : object
            The partially fitted estimator.

        Notes
        -----
        The design matrix is incrementally updated by persisting the gram matrix
        (D.T @ D) and the moment vector (D.T @ y) for each batch, then adding
        the gram and moment contributions of each new batch.

        For batches D = [D1, D2, ..., Dk].T::

            [D1; D2; ...; Dk].T @ [D1; D2; ...; Dk] = sum_k(Dk.T @ Dk)
            [D1; D2; ...; Dk].T @ [y1; y2; ...; yk] = sum_k(Dk.T @ yk)

        This allows incremental accumulation.

        The way we accumulate information across batches prevents us from using
        ``Ridge()`` as we do in full fit. In ``partial_fit()`` we never have access to
        the full design matrix, we only persist the aggregate quantities
        ``D.T @ D``, ``D.T @ y``. Instead, for the regularized case, we directly solve
        the system using ``scipy.linalg.solve(A + reg_mat, B)``, which is mathematically
        equivalent to scikit-learn's ``Ridge(solver='cholesky')``.

        One other difference between full fit and partial fit arises in the direct
        solve path. Because ``pinv()`` is acting on ``D.T @ D`` as opposed to just D,
        the condition number is squared. This may require a lower ``rtol`` to maintain
        numerical stability and avoid loss of information.
        """
        # shape: (n_samples, n_features)
        xp, _, device_ = get_namespace_and_device(X)

        first_call = not hasattr(self, "n_features_in_")
        if first_call:
            dtype = [xp.float64, xp.float32]
        else:
            xp_W, _, device_W = get_namespace_and_device(self.W_[0])
            if (xp != xp_W) or (device_ != device_W):
                raise ValueError(
                    f"All batches must use same device and namespace , "
                    f"{(xp_W, device_W)!r}, but got {(xp, device_)!r}."
                )
            dtype = self.W_[0].dtype

        X, y = validate_data(
            self,
            X,
            y,
            dtype=dtype,
            reset=first_call,
            multi_output=True,
            y_numeric=True,
        )

        y = move_to(y, xp=xp, device=device_)
        y = xp.astype(y, X.dtype, copy=False)
        if y.ndim == 1:
            y = xp.expand_dims(y, axis=-1)

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
        check_same_namespace(X, self, attribute="coef_", method="predict")
        xp, _ = get_namespace(X)

        X = validate_data(self, X, reset=False, dtype=[xp.float64, xp.float32])
        out = super()._predict(X)
        return out

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.array_api_support = True
        return tags
