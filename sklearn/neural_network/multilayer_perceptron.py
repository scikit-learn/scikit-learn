"""Multi-layer Perceptron
"""

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# Licence: BSD 3 clause

import numpy as np

from abc import ABCMeta, abstractmethod
from scipy.optimize import fmin_l_bfgs_b

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..externals import six
from ..preprocessing import LabelBinarizer
from ..utils import gen_even_slices
from ..utils import shuffle
from ..utils import atleast2d_or_csr, check_random_state, column_or_1d
from ..utils.extmath import logistic_sigmoid, safe_sparse_dot


def _identity(X):
    """returns the same input array."""
    return X


def _d_logistic(sigm_X):
    """Implements the derivative of the logistic function.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
    """
    return sigm_X * (1 - sigm_X)


def _softmax(Z):
    """Implements the K-way softmax, (exp(Z).T / exp(Z).sum(axis=1)).T

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
    """
    exp_Z = np.exp(Z.T - Z.max(axis=1)).T
    return (exp_Z.T / exp_Z.sum(axis=1)).T


def _tanh(X):
    """Implements the hyperbolic tan function

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
    """
    return np.tanh(X, X)


def _d_tanh(Z):
    """Implements the derivative of the hyperbolic tan function

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
    """
    Z *= Z
    Z *= -1
    Z += 1
    return Z


def _squared_loss(Y, Z):
    """Implements the square loss for regression."""
    return np.sum((Y - Z) ** 2) / (2 * len(Y))


def _log_loss(Y, Z):
    """Implements Logistic loss for binary class.

    Max/Min clipping is enabled to prevent
    invalid zero value in log computation.
    """
    Z = np.clip(Z, 0.00000001, 0.99999999)
    return -np.sum(Y * np.log(Z) +
                  (1 - Y) * np.log(1 - Z)) / Z.shape[0]


class BaseMultilayerPerceptron(six.with_metaclass(ABCMeta, BaseEstimator)):

    """Base class for MLP classification and regression.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """
    activation_functions = {
        'tanh': _tanh,
        'logistic': logistic_sigmoid,
        'softmax': _softmax
    }
    derivative_functions = {
        'tanh': _d_tanh,
        'logistic': _d_logistic
    }
    loss_functions = {
        'squared_loss': _squared_loss,
        'log_loss': _log_loss,
    }

    @abstractmethod
    def __init__(
        self, n_hidden, activation, algorithm,
            alpha, batch_size, learning_rate, eta0, power_t,
            max_iter, shuffle, random_state, tol, verbose, warm_start):
        self.activation = activation
        self.algorithm = algorithm
        self.alpha = alpha
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.eta0 = eta0
        self.power_t = power_t
        self.max_iter = max_iter
        self.n_hidden = n_hidden
        self.shuffle = shuffle
        self.random_state = random_state
        self.tol = tol
        self.verbose = verbose
        self.warm_start = warm_start

        self.coef_hidden_ = None
        self.t_ = None
        self.eta_ = None

    def _pack(self, W1, W2, b1, b2):
        """Pack the coefficients and intercepts from theta."""
        return np.hstack((W1.ravel(), W2.ravel(),
                          b1.ravel(), b2.ravel()))

    def _unpack(self, theta):
        """Extracts the coefficients and intercepts from theta.

        Parameters
        ----------
        theta : array-like, shape (size(W1) * size(W2) * size(b1) * size(b2))
        A vector comprising the  flattened weights : "W1, W2, b1, b2"
        """
        N = self.n_hidden * self.n_features
        N2 = self.n_outputs * self.n_hidden

        self.coef_hidden_ = np.reshape(
            theta[:N], (self.n_features, self.n_hidden))
        self.coef_output_ = np.reshape(
            theta[N:N2 + N], (self.n_hidden, self.n_outputs))
        self.intercept_hidden_ = theta[N2 + N:N2 + N + self.n_hidden]
        self.intercept_output_ = theta[N2 + N + self.n_hidden:]

    def _validate_params(self):
        """Validate input params. """
        if not isinstance(self.shuffle, bool):
            raise ValueError("shuffle must be either True or False")
        if self.max_iter <= 0:
            raise ValueError("max_iter must be > zero")
        if self.alpha < 0.0:
            raise ValueError("alpha must be >= 0")
        if self.learning_rate in ("constant", "invscaling"):
            if self.eta0 <= 0.0:
                raise ValueError("eta0 must be > 0")

        # raises ValueError if not registered
        if self.activation not in self.activation_functions:
            raise ValueError("The activation %s"
                             " is not supported. " % self.activation)
        if self.learning_rate not in ["constant", "invscaling"]:
            raise ValueError("learning rate %s "
                             " is not supported. " % self.learning_rate)
        if self.algorithm not in ["sgd", "l-bfgs"]:
            raise ValueError("The algorithm %s"
                             " is not supported. " % self.algorithm)

    def _init_fit(self):
        """Initialize weight and bias parameters."""
        rng = check_random_state(self.random_state)

        self.coef_hidden_ = rng.uniform(
            -1, 1, (self.n_features, self.n_hidden))
        self.coef_output_ = rng.uniform(-1, 1, (self.n_hidden, self.n_outputs))
        self.intercept_hidden_ = rng.uniform(-1, 1, self.n_hidden)
        self.intercept_output_ = rng.uniform(-1, 1, self.n_outputs)

    def _init_param(self):
        """Sets the activation, derivative, loss and output functions."""
        self.activation_func = self.activation_functions[self.activation]
        self.derivative_func = self.derivative_functions[self.activation]

        # output for regression
        if self.classes_ is None:
            self.output_func = _identity
        # output for multi class
        elif len(self.classes_) > 2 and self.multi_label is False:
            self.output_func = _softmax
        # output for binary class and multi-label
        else:
            self.output_func = logistic_sigmoid

    def _init_t_eta_(self):
        """Initialize iteration counter attr ``t_``"""
        self.t_ = 1.0
        self.eta_ = self.eta0

    def _preallocate_memory(self, size):
        """ preallocate memory"""
        a_hidden = np.empty((size, self.n_hidden))
        a_output = np.empty((size, self.n_outputs))
        delta_o = np.empty((size, self.n_outputs))

        return a_hidden, a_output, delta_o

    def fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)
             Subset of the target values.

        Returns
        -------
        self
        """
        X = atleast2d_or_csr(X)

        self._validate_params()

        n_samples, self.n_features = X.shape
        self.n_outputs = y.shape[1]

        if not self.warm_start:
            self._init_t_eta_()
            self._init_fit()
            self._init_param()
        else:
            if self.t_ is None or self.coef_hidden_ is None:
                self._init_t_eta_()
                self._init_fit()
                self._init_param()

        if self.shuffle:
            X, y = shuffle(X, y, random_state=self.random_state)

        # l-bfgs does not use mini-batches
        if self.algorithm == 'l-bfgs':
            batch_size = n_samples
        else:
            batch_size = np.clip(self.batch_size, 0, n_samples)
            n_batches = n_samples / batch_size
            batch_slices = list(
                gen_even_slices(
                    n_batches * batch_size,
                    n_batches))

        # preallocate memory
        a_hidden, a_output, delta_o = self._preallocate_memory(
            batch_size)

        if self.algorithm == 'sgd':
            prev_cost = np.inf

            for i in xrange(self.max_iter):
                for batch_slice in batch_slices:
                    cost = self._backprop_sgd(
                        X[batch_slice], y[batch_slice],
                        batch_size, a_hidden,
                        a_output, delta_o)

                if self.verbose:
                    print("Iteration %d, cost = %.2f"
                          % (i, cost))
                if abs(cost - prev_cost) < self.tol:
                    break
                prev_cost = cost
                self.t_ += 1

        elif 'l-bfgs':
            self._backprop_lbfgs(
                X, y, n_samples, a_hidden,
                a_output, delta_o)

        return self

    def _backprop(self, X, y, n_samples, a_hidden, a_output, delta_o):
        """Computes the MLP cost  function
        and its corresponding derivatives with respect to the
        different parameters given in the initialization.

        Parameters
        ----------
        theta : array-like, shape (size(W1) * size(W2) * size(b1) * size(b2))
                    A vector comprising the  flattened weights :
                    "W1, W2, b1, b2"

        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)
            Subset of the target values.

        n_samples : int
            Number of samples

        Returns
        -------
        cost : float
        grad : array-like, shape (size(W1) * size(W2) * size(b1) * size(b2))
        """
        # Forward propagate
        a_hidden[:] = self.activation_func(safe_sparse_dot(X,
                                                           self.coef_hidden_) +
                                           self.intercept_hidden_)

        a_output[:] = self.output_func(safe_sparse_dot(a_hidden,
                                                       self.coef_output_) +
                                       self.intercept_output_)

        # get cost
        cost = self.loss_functions[self.loss](y, a_output)
        # add regularization term to cost
        cost += (0.5 * self.alpha) * (np.sum(self.coef_hidden_ ** 2) +
                                      np.sum(self.coef_output_ ** 2)) \
            / n_samples

        # backward propagate
        diff = y - a_output
        delta_o[:] = -diff
        delta_h = np.dot(delta_o, self.coef_output_.T) *\
            self.derivative_func(a_hidden)

        # get regularized gradient
        W1grad = (safe_sparse_dot(X.T,
                                  delta_h) + (self.alpha *
                                              self.coef_hidden_)) / n_samples
        W2grad = (safe_sparse_dot(a_hidden.T,
                                  delta_o) + (self.alpha *
                                              self.coef_output_)) / n_samples
        b1grad = np.mean(delta_h, 0)
        b2grad = np.mean(delta_o, 0)

        return cost, W1grad, W2grad, b1grad, b2grad

    def _backprop_sgd(
            self, X, y, n_samples, a_hidden, a_output, delta_o):
        """Updates the weights using the computed gradients.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)
            Subset of the target values.

        n_samples : int
            Number of samples.

        """
        cost, W1grad, W2grad, b1grad, b2grad = self._backprop(
            X, y, n_samples, a_hidden, a_output, delta_o)

        # Update weights
        self.coef_hidden_ -= (self.eta_ * W1grad)
        self.coef_output_ -= (self.eta_ * W2grad)
        self.intercept_hidden_ -= (self.eta_ * b1grad)
        self.intercept_output_ -= (self.eta_ * b2grad)

        if self.learning_rate == 'invscaling':
            self.eta_ = self.eta0 / pow(self.t_, self.power_t)

        return cost

    def _backprop_lbfgs(self, X, y, n_samples,
                        a_hidden, a_output, delta_o):
        """Applies the quasi-Newton optimization methods that uses a l_BFGS
        to train the weights.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)
            Subset of the target values.

        n_samples : int
            Number of samples.

        """
        initial_theta = self._pack(
            self.coef_hidden_,
            self.coef_output_,
            self.intercept_hidden_,
            self.intercept_output_)

        if self.verbose is True or self.verbose >= 1:
            iprint = 1
        else:
            iprint = -1

        optTheta, _, _ = fmin_l_bfgs_b(
            func=self._cost_grad,
            x0=initial_theta,
            maxfun=self.max_iter,
            iprint=iprint,
            pgtol=self.tol,
            args=(
                X,
                y,
                n_samples,
                a_hidden, a_output, delta_o))

        self._unpack(optTheta)

    def _cost_grad(self, theta, X, y, n_samples, a_hidden, a_output, delta_o):
        """Computes the MLP cost  function and its
        corresponding derivatives with respect to the
        different parameters given in the initialization.

        Parameters
        ----------
        theta : array-like, shape (size(W1) * size(W2) * size(b1) * size(b2))
                    A vector comprising the  flattened weights :
                    "W1, W2, b1, b2"

        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)
            Subset of the target values.

        n_samples : int
            Number of samples.

        Returns
        -------
        cost : float
        grad : array-like, shape (size(W1) * size(W2) * size(b1) * size(b2))

        """
        self._unpack(theta)

        cost, W1grad, W2grad, b1grad, b2grad = self._backprop(
            X, y, n_samples, a_hidden, a_output, delta_o)
        grad = self._pack(W1grad, W2grad, b1grad, b2grad)

        return cost, grad

    def partial_fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Subset of training data.

        y : numpy array of shape (n_samples)
            Subset of target values.

        Returns
        -------
        self : returns an instance of self.
        """
        X = atleast2d_or_csr(X)

        self.n_outputs = y.shape[1]

        n_samples, self.n_features = X.shape
        self._validate_params()

        if self.coef_hidden_ is None:
            self._init_fit()
            self._init_param()

        if self.t_ is None or self.eta_ is None:
            self._init_t_eta_()

        a_hidden, a_output, delta_o = self._preallocate_memory(n_samples)

        cost = self._backprop_sgd(X, y, n_samples, a_hidden, a_output, delta_o)
        if self.verbose:
            print("Iteration %d, cost = %.2f" % (self.t_, cost))
        self.t_ += 1

        return self

    def decision_function(self, X):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)

        Returns
        -------
        array, shape (n_samples)
        Predicted target values per element in X.
        """
        X = atleast2d_or_csr(X)

        a_hidden = self.activation_func(safe_sparse_dot(X, self.coef_hidden_) +
                                        self.intercept_hidden_)
        output = safe_sparse_dot(a_hidden, self.coef_output_) +\
            self.intercept_output_
        if output.shape[1] == 1:
            output = output.ravel()

        return output


class MultilayerPerceptronClassifier(BaseMultilayerPerceptron,
                                     ClassifierMixin):

    """Multi-layer Perceptron classifier.

    Under a loss function, the algorithm trains either by l-bfgs or gradient 
    descent. The training is done iteratively, in that at each time step the 
    partial derivatives of the loss function with respect to the model 
    parameters are used to update the parameters.  

    It has a regularizer as a penalty term added to the loss function that 
    shrinks model parameters towards zero.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    Please note that this implementation uses one hidden layer only.

    Parameters
    ----------
    n_hidden : int, default 100
        Number of units in the hidden layer.

    activation : {'logistic', 'tanh'}, default 'logistic'
        Activation function for the hidden layer.

         - 'logistic' for 1 / (1 + exp(x)).

         - 'tanh' for the hyperbolic tangent.

    algorithm : {'l-bfgs', 'sgd'}, default 'l-bfgs'
        The algorithm for weight optimization.  Defaults to 'l-bfgs'

        - 'l-bfgs' is an optimization algorithm in the family of 
          quasi-Newton methods.

        - 'sgd' refers to stochastic gradient descent.

    alpha : float, optional, default 0.00001
        L2 penalty (regularization term) parameter.

    batch_size : int, optional, default 200
        Size of minibatches in SGD optimizer.
        If you select the algorithm as 'l-bfgs',
        then the classifier will not use minibatches.

    learning_rate : {'constant', 'invscaling'}, default 'constant'
        Base learning rate for weight updates.

        -'constant', as it stands,  keeps the learning rate 'eta' constant
          throughout training. eta = eta0

        -'invscaling' gradually decreases the learning rate 'eta' at each
          time step 't' using an inverse scaling exponent of 'power_t'.
          eta = eta0 / pow(t, power_t)

    max_iter : int, optional, default 200
        Maximum number of iterations. The algorithm
        iterates until convergence (determined by 'tol') or
        this number of iterations.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    shuffle : bool, optional, default False
        Whether to shuffle samples in each iteration before extracting
        minibatches.

    tol : float, optional, default 1e-5
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached and the algorithm exits.

    eta0 : double, optional, default 0.5
        The initial learning rate used. It controls the step-size
        in updating the weights.

    power_t : double, optional, default 0.25
        The exponent for inverse scaling learning rate.
        It is used in updating eta0 when the learning_rate
        is set to 'invscaling'.

    verbose : bool, optional, default False
        Whether to print progress messages to stdout.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.
    """

    def __init__(
            self, n_hidden=100, activation="logistic",
            algorithm='l-bfgs', alpha=0.00001,
            batch_size=200, learning_rate="constant", eta0=0.5,
            power_t=0.25, max_iter=200, shuffle=False,
            random_state=None, tol=1e-5,
            verbose=False, warm_start=False):

        sup = super(MultilayerPerceptronClassifier, self)
        sup.__init__(n_hidden, activation,
                     algorithm, alpha,
                     batch_size, learning_rate,
                     eta0, power_t,
                     max_iter, shuffle,
                     random_state,
                     tol, verbose,
                     warm_start)

        self.loss = 'log_loss'
        self.classes_ = None

    def fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)

        Returns
        -------
        self
        """
        y = column_or_1d(y, warn=True)

        # needs a better way to check multi-label instances
        if isinstance(np.reshape(y, (-1, 1))[0][0], list):
            self.multi_label = True
        else:
            self.multi_label = False

        self.classes_ = np.unique(y)
        self._lbin = LabelBinarizer()
        y = self._lbin.fit_transform(y)

        super(MultilayerPerceptronClassifier, self).fit(X, y)

        return self

    def partial_fit(self, X, y, classes=None):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        classes : array, shape (n_classes)
            Classes across all calls to partial_fit.
            Can be obtained by via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.
            Note that y doesn't need to contain all labels in `classes`.

        y : numpy array of shape (n_samples)
             Subset of the target values.

        Returns
        -------
        self
        """
        if self.algorithm != 'sgd':
            raise ValueError("only SGD algorithm"
                             " supports partial fit")

        if self.classes_ is None and classes is None:
            raise ValueError("classes must be passed on the first call "
                             "to partial_fit.")
        elif self.classes_ is not None and classes is not None:
            if np.any(self.classes_ != np.unique(classes)):
                raise ValueError("`classes` is not the same as on last call "
                                 "to partial_fit.")
        elif classes is not None:
            self.classes_ = classes

        if not hasattr(self, '_lbin'):
            self._lbin = LabelBinarizer()
            self._lbin._classes = classes

        y = column_or_1d(y, warn=True)

        # needs a better way to check multi-label instances
        if isinstance(np.reshape(y, (-1, 1))[0][0], list):
            self.multi_label = True
        else:
            self.multi_label = False

        y = self._lbin.fit_transform(y)
        super(MultilayerPerceptronClassifier, self).partial_fit(X, y)

        return self

    def predict(self, X):
        """Predict using the multi-layer perceptron model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)

        Returns
        -------
        array, shape (n_samples)
            Predicted target values per element in X.
        """
        X = atleast2d_or_csr(X)
        scores = self.decision_function(X)

        if len(scores.shape) == 1 or self.multi_label is True:
            scores = logistic_sigmoid(scores)
            results = (scores > 0.5).astype(np.int)

            if self.multi_label:
                return self._lbin.inverse_transform(results)

        else:
            scores = _softmax(scores)
            results = scores.argmax(axis=1)

        return self.classes_[results]

    def predict_log_proba(self, X):
        """Returns the log of probability estimates.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        T : array-like, shape (n_samples, n_outputs)
            Returns the log-probability of the sample for each class in the
            model, where classes are ordered as they are in
            `self.classes_`. Equivalent to log(predict_proba(X))
        """
        return np.log(self.predict_proba(X))

    def predict_proba(self, X):
        """Probability estimates.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)

        Returns
        -------
        array, shape (n_samples, n_outputs)
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in `self.classes_`.
        """
        scores = self.decision_function(X)

        if len(scores.shape) == 1:
            scores = logistic_sigmoid(scores)
            return np.vstack([1 - scores, scores]).T
        else:
            return _softmax(scores)


class MultilayerPerceptronRegressor(BaseMultilayerPerceptron, RegressorMixin):

    """Multi-layer Perceptron regressor.

    Under a loss function, the algorithm trains either by l-bfgs or gradient 
    descent. The training is done iteratively, in that at each time step the 
    partial derivatives of the loss function with respect to the model 
    parameters are used to update the parameters.  

    It has a regularizer as a penalty term added to the loss function that 
    shrinks model parameters towards zero.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    Please note that this implementation uses one hidden layer only.

    Parameters
    ----------
    n_hidden : int, default 100
        Number of units in the hidden layer.

    activation : {'logistic', 'tanh'}, default 'logistic'
        Activation function for the hidden layer.

         - 'logistic' for 1 / (1 + exp(x)).

         - 'tanh' for the hyperbolic tangent.

    algorithm : {'l-bfgs', 'sgd'}, default 'l-bfgs'
        The algorithm for weight optimization.  Defaults to 'l-bfgs'

        - 'l-bfgs' is an optimization algorithm in the family of 
          quasi-Newton methods.

        - 'sgd' refers to stochastic gradient descent.

    alpha : float, optional, default 0.00001
        L2 penalty (regularization term) parameter.

    batch_size : int, optional, default 200
        Size of minibatches in SGD optimizer.
        If you select the algorithm as 'l-bfgs',
        then the classifier will not use minibatches.

    learning_rate : {'constant', 'invscaling'}, default 'constant'
        Base learning rate for weight updates.

        -'constant', as it stands,  keeps the learning rate 'eta' constant
          throughout training. eta = eta0

        -'invscaling' gradually decreases the learning rate 'eta' at each
          time step 't' using an inverse scaling exponent of 'power_t'.
          eta = eta0 / pow(t, power_t)

    max_iter : int, optional, default 200
        Maximum number of iterations. The algorithm
        iterates until convergence (determined by 'tol') or
        this number of iterations.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    shuffle : bool, optional, default False
        Whether to shuffle samples in each iteration before extracting
        minibatches.

    tol : float, optional, default 1e-5
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached and the algorithm exits.

    eta0 : double, optional, default 0.1
        The initial learning rate used. It controls the step-size
        in updating the weights.

    power_t : double, optional, default 0.25
        The exponent for inverse scaling learning rate.
        It is used in updating eta0 when the learning_rate
        is set to 'invscaling'.

    verbose : bool, optional, default False
        Whether to print progress messages to stdout.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.
    """

    def __init__(
            self, n_hidden=100, activation="logistic",
            algorithm='l-bfgs', alpha=0.00001,
            batch_size=200, learning_rate="constant", eta0=0.1,
            power_t=0.25, max_iter=100, shuffle=False,
            random_state=None, tol=1e-5,
            verbose=False, warm_start=False):

        sup = super(MultilayerPerceptronRegressor, self)
        sup.__init__(n_hidden, activation,
                     algorithm, alpha,
                     batch_size, learning_rate,
                     eta0, power_t,
                     max_iter, shuffle,
                     random_state,
                     tol, verbose,
                     warm_start)

        self.loss = 'squared_loss'
        self.classes_ = None

    def fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)
            Subset of the target values.

        Returns
        -------
        self
        """
        y = np.atleast_1d(y)

        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))

        super(MultilayerPerceptronRegressor, self).fit(X, y)
        return self

    def partial_fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape (n_samples)
            Subset of the target values.

        Returns
        -------
        self
        """
        y = np.atleast_1d(y)

        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))

        super(MultilayerPerceptronRegressor, self).partial_fit(X, y)
        return self

    def predict(self, X):
        """Predict using the multi-layer perceptron model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)

        Returns
        -------
        array, shape (n_samples)
            Predicted target values per element in X.
        """
        X = atleast2d_or_csr(X)

        return self.decision_function(X)
