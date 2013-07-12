"""Sparse Autoencoder
"""

# Author: Issam Laradji <issam.laradji@gmail.com>

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from scipy.linalg import norm
from itertools import cycle, izip
from sklearn.utils import atleast2d_or_csr, check_random_state
from sklearn.utils import gen_even_slices
from sklearn.utils import shuffle
from sklearn.base import BaseEstimator, TransformerMixin


def binary_KL_divergence(p, p_hat):
    """
    Computes the a real, KL divergence of two binomial distributions with
    probabilities p  and p_hat respectively.
    """
    return (p * np.log(p / p_hat)) + ((1 - p) * np.log((1 - p) / (1 - p_hat)))


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


class Autoencoder(BaseEstimator, TransformerMixin):

    """
    Sparse Autoencoder (SAE)

    A Sparse Autoencoder with one hidden layer.
    Parameters
    ----------
    n_hidden : int
        Number of hidden neurons
    activation: string, optional
        Activation function for the hidden layer; either "logistic" for
        1 / (1 + exp(x)), or "tanh" for the hyperbolic tangent.
    algorithm : string, optional
        Optimization function for training the weights; could be "l-bfgs-b", "cg",
        "newton-cg", or "bfgs"
    learning_rate : float, optional
        Learning rate to use during learning. It is *highly* recommended
        to tune this hyper-parameter. Possible values are 10**[0., -3.].
    beta : float, optional
        Weight of sparsity penalty term
    sparsity_param : float, optional
        Desired average activation of the hidden units
    batch_size : int, optional
        Number of examples per minibatch.
    max_iter : int, optional
        Number of iterations/sweeps over the training dataset to perform
        during training.
    tol : float, optional
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached.
    verbose: bool, optional
        When True (False by default) the method outputs the progress
        of learning after each iteration.
    random_state : integer or numpy.RandomState, optional
        A random number generator instance to define the state of the
        random permutations generator. If an integer is given, it fixes the
        seed. Defaults to the global numpy random number generator.

    Attributes
    ----------
    self.coef_hidden_ : array-like, shape (n_hidden, n_features)
        Weight matrix, where n_features in the number of visible
        units and n_hidden is the number of hidden units.
    self.coef_output_  : array-like, shape (n_features, n_hidden)
        Weight matrix, where n_features in the number of visible
        units and n_hidden is the number of hidden units.
    intercept_hidden_ : array-like, shape (n_hidden,), optional
        Biases of the hidden units
    intercept_visible_ : array-like, shape (n_features,), optional
        Biases of the visible units

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.neural_network import SAE
    >>> X = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]])
    >>> model = SAE(n_hidden=10)
    >>> model.fit(X)
    Autoencoder(activation_func='logistic', alpha=0.0001, batch_size=1000, beta=3,
  learning_rate=0.0001, max_iter=20, n_hidden=10,
  algorithm='l-bfgs', random_state=None, sparsity_param=0.01,
  tol=1e-05, verbose=False)

    References
    ----------

    [1] Ngiam, Jiquan, et al. "On optimization methods for deep learning."
        Proceedings of the 28th International Conference on Machine Learning (ICML-11). 2011.
        http://ai.stanford.edu/~quocle/LeNgiCoaLahProNg11.pdf
    """
    activation_functions = {
            'tanh': tanh,
            'logistic': logistic
        }
    derivative_functions = {
            'tanh': d_tanh,
            'logistic': d_logistic
        }
    def __init__(
        self, n_hidden=25, activation_func='logistic', algorithm='l-bfgs',
        learning_rate=0.3, alpha=3e-3, beta=3, sparsity_param=0.1,
            batch_size=100, shuffle_data=False, max_iter=20, tol=1e-5, verbose=False, random_state=None):
        self.activation_func = activation_func
        self.algorithm = algorithm
        self.n_hidden = n_hidden
        self.alpha = alpha
        self.learning_rate = learning_rate
        self.beta = beta
        self.sparsity_param = sparsity_param
        self.batch_size = batch_size
        self.shuffle_data = shuffle_data
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose
        self.random_state = random_state

    def _init_fit(self, n_features):
        """
        Initialize weight and bias parameters

        Parameters
        ----------
        n_features: int
            Number of features (visible nodes).

        Returns
        -------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
        """
        rng = check_random_state(self.random_state)
        self.coef_hidden_ = rng.uniform(-1, 1, (n_features, self.n_hidden))
        self.coef_output_ = rng.uniform(-1, 1, (self.n_hidden, n_features))
        self.intercept_hidden_ = rng.uniform(-1, 1, self.n_hidden)
        self.intercept_output_ = rng.uniform(-1, 1, n_features)

    def _init_param(self):
        """
        Sets the activation, derivative and the output functions
        """
        self.activation = self.activation_functions[self.activation_func]
        self.derivative = self.derivative_functions[self.activation_func]
        
    def _unpack(self, theta, n_features):
        """
        Extract the coefficients and intercepts (W1,W2,b1,b2) from theta

        Parameters
        ----------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
          Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
        n_features: int
          Number of features (visible nodes).
        """
        N = self.n_hidden * n_features
        self.coef_hidden_ = np.reshape(theta[:N],
                                      (n_features, self.n_hidden))
        self.coef_output_ = np.reshape(theta[N:2 * N],
                                      (self.n_hidden, n_features))
        self.intercept_hidden_ = theta[2 * N:2 * N + self.n_hidden]
        self.intercept_output_ = theta[2 * N + self.n_hidden:]

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

    def transform(self, X):
        """
        Computes the extracted features.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        Returns
        -------
        h: array-like, shape (n_samples, n_components)
        """
        return self.activation(np.dot(X, self.coef_hidden_) + self.intercept_hidden_)

    def fit_transform(self, X):
        """
        Fit the model to the data X and transform it.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        """
        self.fit(X)
        return self.transform(X)

    def fit(self, X):
        """
        Fit the model to the data X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self
        """
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        n_samples, n_features = X.shape
        self._init_fit(n_features)
        self._init_param()
        if self.shuffle_data:
            X, y = shuffle(X, y, random_state=self.random_state)
        # generate batch slices
        self.batch_size = np.clip(self.batch_size, 0, n_samples)
        n_batches = n_samples / self.batch_size
        batch_slices = list(
            gen_even_slices(
                n_batches *
                self.batch_size,
                n_batches))
        #l-bfgs does not work well with minibatches
        if self.algorithm == 'l-bfgs':
            self.batch_size = n_samples
        # preallocate memory
        a_hidden = np.empty((self.batch_size, self.n_hidden))
        a_output = np.empty((self.batch_size, n_features))
        delta_o = np.empty((self.batch_size, n_features))
        if self.algorithm == 'sgd':
            for i in xrange(self.max_iter):
                    for batch_slice in batch_slices:
                        cost = self.backprop_sgd(
                            X[batch_slice],
                            n_features, self.batch_size,
                            delta_o, a_hidden, a_output)
                    if self.verbose:
                        print("Iteration %d, cost = %.2f"
                              % (i, cost))
        elif self.algorithm == 'l-bfgs':
            self._backprop_lbfgs(
                X, n_features,
                a_hidden, a_output, 
                delta_o, n_samples)
        return self

    def backprop(self, X, n_features, n_samples,
                  delta_o, a_hidden, a_output):
        """
        Computes the sparse autoencoder cost  function ``Jsparse(W,b)``
        and the corresponding derivatives of Jsparse with respect to the
        different parameters given in the initialization [1]

        Parameters
        ----------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))
          Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
        X: array-like, shape (n_samples, n_features)
          Training data, where n_samples in the number of samples
          and n_features is the number of features.
        n_features: int
          Number of features (visible nodes).
        n_samples: int
          Number of samples

       Returns
       -------
       cost: float
       grad: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))

       References
       -------
       [1] http://ufldl.stanford.edu/wiki/index.php/Autoencoders_and_Sparsity
        """
        # Forward propagate
        a_hidden[:] = self.activation(np.dot(X, self.coef_hidden_)
                                      + self.intercept_hidden_)
        a_output[:] = self.activation(np.dot(a_hidden, self.coef_output_)
                                      + self.intercept_output_)
        # Get average activation of hidden neurons
        sparsity_param_hat = np.sum(a_hidden, 0) / n_samples
        sparsity_delta  = self.beta * \
            ((1 - self.sparsity_param) / (1 - sparsity_param_hat)
             - self.sparsity_param / sparsity_param_hat)
        # Backward propagate
        diff = X - a_output
        delta_o[:] = -diff * self.derivative(a_output)
        delta_h = (
            (np.dot(delta_o, self.coef_output_.T) +
             sparsity_delta)) *\
            self.derivative(a_hidden)
        # Get cost 
        cost = np.sum(diff ** 2) / (2 * n_samples)
        # Add regularization term to cost 
        cost += (0.5 * self.alpha) * (
            np.sum(self.coef_hidden_ ** 2) + np.sum(
                self.coef_output_ ** 2))
        # Add sparsity term to the cost
        cost += self.beta * np.sum(
            binary_KL_divergence(
                self.sparsity_param,
                sparsity_param_hat))
        #Get gradients
        W1grad = np.dot(X.T, delta_h) / n_samples 
        W2grad = np.dot(a_hidden.T, delta_o) / n_samples
        b1grad = np.sum(delta_h, 0) / n_samples
        b2grad = np.sum(delta_o, 0) / n_samples
        # Add regularization term to gradients 
        W1grad += self.alpha * self.coef_hidden_
        W2grad += self.alpha * self.coef_output_
        return cost, W1grad, W2grad, b1grad, b2grad

    def backprop_sgd(
            self, X, n_features, n_samples, delta_o, a_hidden, a_output):
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
            X, n_features, n_samples, delta_o, a_hidden, a_output)
        # Update weights
        self.coef_hidden_ -= (self.learning_rate * W1grad)
        self.coef_output_ -= (self.learning_rate * W2grad)
        self.intercept_hidden_ -= (self.learning_rate * b1grad)
        self.intercept_output_ -= (self.learning_rate * b2grad)
        # TODO: dynamically update learning rate
        return cost
        
    def _backprop_lbfgs(
            self, X, n_features, a_hidden, a_output, delta_o, n_samples):
        """
        Applies the one of the optimization methods (l-bfgs-b, bfgs, newton-cg, cg)
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
                n_features,
                n_samples,
                delta_o,
                a_hidden,
                a_output))
        self._unpack(optTheta, n_features)

    def _cost_grad(self, theta, X, n_features,
                   n_samples, delta_o, a_hidden, a_output):
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
        self._unpack(theta, n_features)
        cost, W1grad, W2grad, b1grad, b2grad = self.backprop(
            X, n_features, n_samples, delta_o, a_hidden, a_output)
        return cost, self._pack(W1grad, W2grad, b1grad, b2grad)
