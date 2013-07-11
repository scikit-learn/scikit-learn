"""Multi-layer perceptron
"""

# Author: Issam Laradji <issam.laradji@gmail.com>
# Credits to: Amueller and Larsmans codes on MLP

import numpy as np
from scipy.linalg import norm

from sklearn.utils import gen_even_slices
from sklearn.utils import shuffle
from scipy.optimize import fmin_l_bfgs_b
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin
from sklearn.preprocessing import LabelBinarizer
from sklearn.utils import atleast2d_or_csr, check_random_state
from sklearn.utils.extmath import logsumexp, safe_sparse_dot
from itertools import cycle, izip


def logistic(x):
    """
    Implements the logistic function.

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    return 1. / (1. + np.exp(np.clip(-x, -30, 30)))


def d_logistic(x):
    """
    Implements the derivative of the logistic function.

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    return x * (1 - x)


def log_softmax(X):
    """
    Computes the logistic K-way softmax, (exp(X).T / exp(X).sum(axis=1)).T,
    in the log domain

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    return (X.T - logsumexp(X, axis=1)).T


def softmax(X):
    """
    Computes the K-way softmax, (exp(X).T / exp(X).sum(axis=1)).T,
    in the log domain

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    exp_X = np.exp(X)
    return (exp_X.T / exp_X.sum(axis=1)).T


def tanh(X):
    """
    Computes the hyperbolic tan function

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    return np.tanh(X, X)


def d_tanh(X):
    """
    Computes the derivative of the hyperbolic tan function

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    X *= -X
    X += 1
    return X


class BaseMLP(BaseEstimator):

    """Multi-layer perceptron (feedforward neural network) classifier.

    Trained with gradient descent under square loss function

    Parameters
    ----------
    n_hidden : int
        Number of units in the hidden layer.
    activation: string, optional
        Activation function for the hidden layer; either "logistic" for
        1 / (1 + exp(x)), or "tanh" for the hyperbolic tangent.
    alpha : float, optional
        L2 penalty (weight decay) parameter.
    batch_size : int, optional
        Size of minibatches in SGD optimizer.
    learning_rate : float, optional
        Base learning rate. This will be scaled by sqrt(n_features) for the
        input-to-hidden weights, and by sqrt(n_hidden) for the hidden-to-output
        weights.
    max_iter : int, optional
        Maximum number of iterations.
    random_state : int or RandomState, optional
        State of or seed for random number generator.
    shuffle : bool, optional
        Whether to shuffle samples in each iteration before extracting
        minibatches.
    tol : float, optional
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached.
    verbose : bool, optional
        Whether to print progress messages to stdout.

    """
    activation_functions = {
            'tanh': tanh,
            'logistic': logistic,
            'softmax': softmax
        }
    derivative_functions = {
            'tanh': d_tanh,
            'logistic': d_logistic
        }
    def __init__(
        self, n_hidden, activation, output_func, loss, algorithm,
            alpha, batch_size, learning_rate, eta0, power_t,
            max_iter, shuffle_data, random_state, tol, verbose):
        self.activation = activation
        self.output_func = output_func
        self.loss = loss
        self.algorithm = algorithm
        self.alpha = alpha
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.eta0 = eta0
        self.power_t = power_t
        self.max_iter = max_iter
        self.n_hidden = n_hidden
        self.shuffle_data = shuffle_data
        self.random_state = random_state
        self.tol = tol
        self.verbose = verbose

    def _pack(self, W1, W2, b1, b2):
        """
        Pack the coefficients and intercepts (W1,W2,b1,b2) from theta

        Parameters
        ----------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
            Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
        n_features: int
            Number of features
        n_classes: int
            Number of target classes
        """
        return np.hstack((W1.ravel(), W2.ravel(),
                          b1.ravel(), b2.ravel()))

    def _unpack(self, theta, n_features, n_classes):
        """
        Extracts the coefficients and intercepts (W1,W2,b1,b2) from theta

        Parameters
        ----------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
            Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
        n_features: int
            Number of features
        n_classes: int
            Number of target classes
        """
        N = self.n_hidden * n_features
        N2 = n_classes * self.n_hidden
        self.coef_hidden_ = np.reshape(theta[:N], (n_features, self.n_hidden))
        self.coef_output_ = np.reshape(
            theta[N:N2 + N], (self.n_hidden, n_classes))
        self.intercept_hidden_ = theta[N2:N2 + self.n_hidden]
        self.intercept_output_ = theta[N2 + N + self.n_hidden:]

    def _init_fit(self, n_features, n_classes):
        """
        Initialize weight and bias parameters

        Parameters
        ----------
        n_features: int
            Number of features

        n_classes: int
            Number of target classes

        """
        n_hidden = self.n_hidden
        rng = check_random_state(self.random_state)
        self.coef_hidden_ = rng.uniform(-1, 1, (n_features, n_hidden))
        self.coef_output_ = rng.uniform(-1, 1, (n_hidden, n_classes))
        self.intercept_hidden_ = rng.uniform(-1, 1, n_hidden)
        self.intercept_output_ = rng.uniform(-1, 1, n_classes)

    def _init_param(self):
        """
        Sets the activation, derivative and the output functions
        
        """
        self.activation_func = self.activation_functions[self.activation]
        self.derivative_func = self.derivative_functions[self.activation]
        self.output_func =  self.activation_functions[self.output_func]
        
    def fit(self, X, y):
        """
        Fit the model to the data X and target y.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : numpy array of shape [n_samples]
            Subset of the target values.

        Returns
        -------
        self
        """
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        self._lbin = LabelBinarizer()
        Y = self._lbin.fit_transform(y)
        if Y.shape[1] == 1:
            self.output_func = self.activation
        n_samples, n_features = X.shape
        n_classes = Y.shape[1]
        self._init_fit(n_features, n_classes)
        self._init_param()
        if self.shuffle_data:
            X, Y = shuffle(X, Y, random_state=self.random_state)
        self.batch_size = np.clip(self.batch_size, 0, n_samples)
        n_batches = n_samples / self.batch_size
        batch_slices = list(
            gen_even_slices(
                n_batches *
                self.batch_size,
                n_batches))
        # preallocate memory
        a_hidden = np.empty((self.batch_size, self.n_hidden))
        a_output = np.empty((self.batch_size, n_classes))
        delta_o = np.empty((self.batch_size, n_classes))
        if self.algorithm is 'sgd':
            eta = self.eta0
            t = 1
            for i in  xrange(self.max_iter):
                for batch_slice in batch_slices:
                    cost, eta = self.backprop_naive(
                        X[batch_slice],
                        Y[batch_slice],
                        self.batch_size,
                        a_hidden,
                        a_output,
                        delta_o,
                        t,
                        eta)
                if self.verbose:
                        print("Iteration %d, cost = %.2f"
                              % (i, cost))
                t += 1
        elif 'l-bfgs':
            for batch_slice in batch_slices:
                self._backprop_optimizer(
                    X[batch_slice], Y[batch_slice], n_features, n_classes, self.batch_size, a_hidden,
                     a_output,
                    delta_o)
        return self

    def _backprop_optimizer(
            self, X, Y, n_features, n_classes, n_samples,
             a_hidden, a_output, delta_o):
        """
        Applies the quasi-Newton optimization methods that uses a l_BFGS
        to train the weights

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Y : numpy array of shape [n_samples]
            Subset of the target values.

        n_features: int
            Number of features

        n_classes: int
            Number of target classes

        n_samples: int
            Number of samples

        """
        initial_theta = self._pack(
            self.coef_hidden_,
            self.coef_output_,
            self.intercept_hidden_,
            self.intercept_output_)
        optTheta, _, _ = fmin_l_bfgs_b(
            func=self._cost_grad,
            x0=initial_theta,
            maxfun=self.max_iter,
            disp=self.verbose,
            args=(
                X,
                Y,
                n_features,
                n_classes,
                n_samples,
                a_hidden, a_output, delta_o))
        self._unpack(optTheta, n_features, n_classes)

    def _cost_grad(self, theta, X, Y, n_features, n_classes,
                   n_samples, a_hidden, a_output, delta_o):
        """
        Computes the MLP cost  function ``J(W,b)``
        and the corresponding derivatives of J(W,b) with respect to the
        different parameters given in the initialization

        Parameters
        ----------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))
            Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        n_features: int
            Number of features
        n_classes: int
            Number of target classes
        n_samples: int
            Number of samples

        Returns
        -------
        cost: float
        grad: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))

        """
        self._unpack(theta, n_features, n_classes)
        cost, W1grad, W2grad, b1grad, b2grad = self.backprop(
            X, Y, n_samples, a_hidden, a_output, delta_o)
        grad = self._pack(W1grad, W2grad, b1grad, b2grad)
        return cost, grad

    def backprop_naive(
            self, X, Y, n_samples, a_hidden, a_output, delta_o, t, eta):
        """
        Updates the weights using the computed gradients

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Y : numpy array of shape [n_samples]
            Subset of the target values.

        n_features: int
            Number of features

        n_classes: int
            Number of target classes

        n_samples: int
            Number of samples

        """
        cost, W1grad, W2grad, b1grad, b2grad = self.backprop(
            X, Y, n_samples, a_hidden, a_output, delta_o)
        # Update weights
        self.coef_hidden_ -= (eta * W1grad)
        self.coef_output_ -= (eta * W2grad)
        self.intercept_hidden_ -= (eta * b1grad)
        self.intercept_output_ -= (eta * b2grad)
        if self.learning_rate == 'optimal':
            eta = 1.0 / (self.alpha * t)
        elif self.learning_rate == 'invscaling':
            eta = self.eta0 / pow(t, self.power_t)
        return cost, eta

    def backprop(self, X, Y, n_samples, a_hidden, a_output, delta_o):
        """
        Computes the MLP cost  function ``J(W,b)``
        and the corresponding derivatives of J(W,b) with respect to the
        different parameters given in the initialization

        Parameters
        ----------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))
            Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        n_features: int
            Number of features
        n_classes: int
            Number of target classes
        n_samples: int
            Number of samples

        Returns
        -------
        cost: float
        grad: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))

        """
        # Forward propagate
        a_hidden[:] = self.activation_func(safe_sparse_dot(X, self.coef_hidden_) +
                                      self.intercept_hidden_)
        a_output[:] = self.output_func(safe_sparse_dot(a_hidden, self.coef_output_) +
                                       self.intercept_output_)
        # Backward propagate
        diff = Y - a_output
        if self.loss == 'squared_loss':
            delta_o[:] = -diff * self.derivative_func(a_output)
            delta_h = np.dot(
                delta_o,
                self.coef_output_.T) * self.derivative_func(a_hidden)
            cost = np.sum(diff**2)/ (2 * n_samples)
        elif self.loss == 'log':
            delta_o[:] = -diff
            delta_h = np.dot(delta_o, self.coef_output_.T) *\
                    self.derivative_func(a_hidden)
            cost = np.sum(
                np.sum(-Y * np.log(a_output) - (1 - Y) * np.log(1 - a_output)))
        # Get regularized gradient
        W1grad = safe_sparse_dot(X.T, delta_h) + \
                (self.alpha * self.coef_hidden_)
        W2grad = safe_sparse_dot(a_hidden.T, delta_o) + \
                (self.alpha * self.coef_output_)
        b1grad = np.mean(delta_h, 0)
        b2grad = np.mean(delta_o, 0)
        # Add regularization term to cost
        cost += (
            0.5 * self.alpha) * (
                np.sum(
                    self.coef_hidden_ ** 2) + np.sum(
                        self.coef_output_ ** 2))
        W1grad /= n_samples
        W2grad /= n_samples
        return cost, W1grad, W2grad, b1grad, b2grad

    def partial_fit(self, X, y, classes):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
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
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        _, n_features = X.shape
        self._init_param()
        if self.classes_ is None and classes is None:
            raise ValueError("classes must be passed on the first call "
                             "to partial_fit.")
        elif classes is not None and self.classes_ is not None:
            if not np.all(self.classes_ == np.unique(classes)):
                raise ValueError("`classes` is not the same as on last call "
                                 "to partial_fit.")
        elif classes is not None:
            self._lbin = LabelBinarizer(classes=classes)
            Y = self._lbin.fit_transform(y)
            self._init_fit(n_features, Y.shape[1])
        else:
            Y = self._lbin.transform(y)

        self.backprop_naive(X, Y, 1)
        return self

    def decision_function(self, X):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples]
           Predicted target values per element in X.
        """
        a_hidden = self.activation_func(safe_sparse_dot(X, self.coef_hidden_) +
                                   self.intercept_hidden_)
        output = self.output_func(safe_sparse_dot(a_hidden, self.coef_output_) +
                                  self.intercept_output_)
        if output.shape[1] == 1:
            output = output.ravel()
        return output

    def predict(self, X):
        """Predict using the MLP model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples]
           Predicted target values per element in X.
        """
        X = atleast2d_or_csr(X)
        scores = self.decision_function(X)
        if len(scores.shape) == 1:
            indices = (scores > 0).astype(np.int)
        else:
            indices = scores.argmax(axis=1)
        return self.classes_[indices]

    def predict_log_proba(self, X):
        """Log of probability estimates.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the log-probability of the sample for each class in the
            model, where classes are ordered as they are in
            `self.classes_`.
        """
        scores = self.decision_function(X)
        if len(scores.shape) == 1:
            return np.log(self.activation_func(scores))
        else:
            return log_softmax(scores)

    def predict_proba(self, X):
        """Probability estimates.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in `self.classes_`.
        """
        prob = self.decision_function(X)
        prob *= -1
        np.exp(prob, prob)
        prob += 1
        np.reciprocal(prob, prob)
        if len(prob.shape) == 1:
            return np.vstack([1 - prob, prob]).T 
        else:
            return softmax(scores)

    @property
    def classes_(self):
        return self._lbin.classes_


class MLPClassifier(BaseMLP, ClassifierMixin):

    def __init__(
        self, n_hidden=100, activation="tanh", output_func='softmax',
        loss='log', algorithm='sgd', alpha=0.00001, batch_size=2000,
        learning_rate="constant", eta0=1, power_t=0.5, max_iter=100,
        shuffle_data=False, random_state=None, tol=1e-5, verbose=False):
        super(
            MLPClassifier, self).__init__(n_hidden, activation, output_func, loss,
                                          algorithm, alpha, batch_size, learning_rate, eta0, 
                                          power_t, max_iter, shuffle_data, random_state, tol, verbose)


class MLPRegressor(BaseMLP, RegressorMixin):

    def __init__(
        self, n_hidden=100, activation="tanh", output_func='tanh',
         loss='squared_loss', algorithm='sgd', alpha=0.00001, 
         batch_size=2000, learning_rate="constant", eta0=0.1, 
         power_t=0.5, max_iter=100, shuffle_data=False, 
         random_state=None, tol=1e-5, verbose=False):
        super(
            MLPClassifier, self).__init__(n_hidden, activation, output_func, loss,
                                          algorithm, alpha, batch_size, learning_rate, eta0, 
                                          power_t, max_iter, shuffle_data, random_state, tol, verbose)
