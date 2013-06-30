"""Mulit-layer perceptron
"""

# Author: Issam Laradji <issam.laradji@gmail.com>
# Credits to: Amueller and Larsmans codes on MLP

import numpy as np
from scipy.linalg import norm
from scipy.optimize import fmin_l_bfgs_b

from ..base import BaseEstimator, ClassifierMixin
from ..preprocessing import LabelBinarizer
from ..utils import atleast2d_or_csr, check_random_state
from ..utils.extmath import logsumexp, safe_sparse_dot


def validate_grad(J, theta, n_slice):
    """
    (For debugging only, it will be removed later)
    Validates the gradient  by comparing between
    analytic and numerical gradient computation
    """
    numgrad = np.zeros(n_slice)
    E = np.eye(np.size(theta))
    epsilon = 1e-5
    for i in range(20):
        dtheta = E[:, i] * epsilon
        numgrad[i] = (J(theta + dtheta) - J(theta - dtheta)) / epsilon / 2.0
    return numgrad

# TODO: Add more activiation functions


def sigmoid(x):
    """
    Implements the sigmoid function.

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    return 1. / (1. + np.exp(np.clip(-x, -30, 30)))


def d_sigmoid(x):
    """
    Implements the derivative of the sigmoid function.

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
    Computes the sigmoid K-way softmax, (exp(X).T / exp(X).sum(axis=1)).T,
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
    Computes the sigmoid K-way softmax, (exp(X).T / exp(X).sum(axis=1)).T,
    in the log domain

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    if X.shape[1] == 1:
        return sigmoid(X)
    else:
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


class MLPClassifier(BaseEstimator, ClassifierMixin):

    """Multi-layer perceptron (feedforward neural network) classifier.

    Trained with gradient descent under square loss function

    Parameters
    ----------
    n_hidden : int
        Number of units in the hidden layer.
    activation: string, optional
        Activation function for the hidden layer; either "sigmoid" for
        1 / (1 + exp(x)), or "tanh" for the hyperbolic tangent.
    _lambda : float, optional
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

    def __init__(
        self, n_hidden, activation="sigmoid", l_bfgs=False,
        _lambda=0, batch_size=100, learning_rate=.1,
        max_iter=100, random_state=None,
            shuffle=True, tol=1e-5, verbose=False):
        # TODO: Move these function initializations somewhere else
        activation_functions = {
            'tanh': tanh,
            'sigmoid': sigmoid
        }
        derivative_functions = {
            'tanh': d_tanh,
            'sigmoid': d_sigmoid
        }
        # TODO: Loss functions
        self.activation = activation_functions[activation]
        self.derivative = derivative_functions[activation]
        self.l_bfgs = l_bfgs
        self._lambda = _lambda
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.max_iter = max_iter
        self.n_hidden = n_hidden
        self.random_state = check_random_state(random_state)
        self.shuffle = shuffle
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
        self.coef_hidden_ = np.reshape(theta[:N], (self.n_hidden, n_features))
        self.coef_output_ = np.reshape(
            theta[N:N2 + N], (n_classes, self.n_hidden))
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

        Returns
        -------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
        """
        n_hidden = self.n_hidden
        rng = self.random_state
        self.coef_hidden_ = rng.uniform(-1, 1, (n_hidden, n_features))
        self.coef_output_ = rng.uniform(-1, 1, (n_classes, n_hidden))
        self.intercept_hidden_ = rng.uniform(-1, 1, n_hidden)
        self.intercept_output_ = rng.uniform(-1, 1, n_classes)

    def benchmark(self, X, Y, n_features, n_classes, n_samples):
        """
        (For debugging only, it will be removed later)
        Showcases the performance of l_bfgs vs. naive backpropagation
        """
        import time
        initial_theta = self._pack(
            self.coef_hidden_,
            self.coef_output_,
            self.intercept_hidden_,
            self.intercept_output_)
        # validate gradient
        print 'Analytically computed gradient'
        print validate_grad(lambda x: self._cost_grad(X, Y, x,
                                                      n_features, n_classes, n_samples)[
                            0],
                            initial_theta, 20)
        print 'Numerically computed gradient'
        print self._cost_grad(X, Y, initial_theta,
                              n_features, n_classes,
                              n_samples)[1][:20]
        begin = time.time()
        self.backprop_l_bfgs(X, Y, n_features, n_classes, n_samples)
        print 'l_bfgs time:', time.time() - begin
        # evaluate iterations
        begin = time.time()
        self._init_fit(n_features, n_classes)
        prev_cost = 1000000
        for i in xrange(self.max_iter):
            cost = self.backprop_naive(X, Y, n_samples)[0]
            if abs(cost - prev_cost) < self.tol:
                print 'convergence reached'
                break
            prev_cost = cost
        print 'naive time: ', time.time() - begin
        print 'naive cost:', cost

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
        rng = check_random_state(self.random_state)
        self._lbin = LabelBinarizer()
        Y = self._lbin.fit_transform(y)
        n_features = X.shape[1]
        n_classes = Y.shape[1]
        n_samples = X.shape[0]
        self._init_fit(n_features, n_classes)
        inds = np.arange(n_samples)
        rng.shuffle(inds)
        n_batches = int(np.ceil(len(inds) / float(self.batch_size)))
        # Transpose improves performance (from 0.5 seconds to 0.05)
        X = X.T
        Y = Y.T
        if self.l_bfgs:
            self.backprop_l_bfgs(X, Y, n_features, n_classes, n_samples)
        else:
            prev_cost = 1000000
            for i in xrange(self.max_iter):
                for minibatch in xrange(n_batches):
                    cost = self.backprop_naive(
                        X[:, inds[minibatch::n_batches]],
                        Y[:, inds[minibatch::n_batches]],
                        n_samples)
                # Check tolerance
                if self.verbose:
                    print("Iteration %d, cost = %.2f"
                          % (i, cost))
                if abs(cost - prev_cost) < self.tol:
                    break
                prev_cost = cost
        return self

    def backprop_l_bfgs(self, X, Y, n_features, n_classes, n_samples):
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
        options = {'maxfun': self.max_iter, 'disp': False, 'pgtol': self.tol}
        optTheta, cost, d = \
            fmin_l_bfgs_b(
                lambda t: self._cost_grad(
                    X, Y, t, n_features,
                    n_classes, n_samples),
                initial_theta, **options)
        self._unpack(optTheta, n_features, n_classes)

    def _cost_grad(self, X, Y, theta, n_features, n_classes, n_samples):
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
        # Unpack theta to extract W1, W2, b1, b2
        self._unpack(theta, n_features, n_classes)
        cost, W1grad, W2grad, b1grad, b2grad = self.backprop(X, Y, n_samples)
        # Convert grad to vector form
        # This is necessary for the fmin optimizer)
        grad = self._pack(W1grad, W2grad, b1grad, b2grad)
        return cost, grad

    def backprop_naive(self, X, Y, n_samples):
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
        cost, W1grad, W2grad, b1grad, b2grad = self.backprop(X, Y, n_samples)
        # Update weights
        self.coef_hidden_ -= (self.learning_rate * W1grad)
        self.coef_output_ -= (self.learning_rate * W2grad)
        self.intercept_hidden_ -= (self.learning_rate * b1grad)
        self.intercept_output_ -= (self.learning_rate * b2grad)
        # TODO: dynamically update learning rate
        return cost

    def backprop(self, X, Y, n_samples):
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
        a2 = self.activation(safe_sparse_dot(self.coef_hidden_, X) +
                             self.intercept_hidden_[:, np.newaxis])
        a3 = self.activation(safe_sparse_dot(self.coef_output_, a2) +
                             self.intercept_output_[:, np.newaxis])
        # Backward propagate
        diff = Y - a3
        delta3 = -diff * self.derivative(a3)
        delta2 = np.dot(self.coef_output_.T, delta3) * self.derivative(a2)
        # Get cost and gradient
        cost = np.trace(safe_sparse_dot(diff, diff.T)) / (2 * n_samples)
        W1grad = safe_sparse_dot(delta2, X.T) / n_samples
        W2grad = safe_sparse_dot(delta3, a2.T) / n_samples
        b1grad = np.mean(delta2, 1)
        b2grad = np.mean(delta3, 1)
        # Add regularization term to cost and gradient
        cost += (self._lambda / 2) * \
                (norm(self.coef_hidden_, 'fro') ** 2 +
                 norm(self.coef_output_, 'fro') ** 2)
        W1grad += (self._lambda * self.coef_hidden_)
        W2grad += (self._lambda * self.coef_output_)
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
        a2 = self.activation(safe_sparse_dot(X, self.coef_hidden_.T) +
                             self.intercept_hidden_)
        output = self.activation(safe_sparse_dot(a2, self.coef_output_.T) +
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
            return np.log(self.activation(scores))
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
        scores = self.decision_function(X)
        if len(scores.shape) == 1:
            return self.activation(scores)
        else:
            return softmax(scores)

    @property
    def classes_(self):
        return self._lbin.classes_
