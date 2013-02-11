# Author: Lars Buitinck <L.J.Buitinck@uva.nl>

import numpy as np

from ..base import BaseEstimator, ClassifierMixin
from ..preprocessing import LabelBinarizer
from ..utils import atleast2d_or_csr, check_random_state
from ..utils.extmath import logsumexp, safe_sparse_dot

from .backprop_sgd import backprop_sgd


def logistic(x):
    return 1. / (1. + np.exp(-x))


def log_softmax(X):
    # Computes the logistic K-way softmax, (exp(X).T / exp(X).sum(axis=1)).T,
    # in the log domain
    return (X.T - logsumexp(X, axis=1)).T


def softmax(X):
    if X.shape[1] == 1:
        return logistic(X)
    else:
        exp_X = np.exp(X)
        return (exp_X.T / exp_X.sum(axis=1)).T


def _tanh(X):
    """Hyperbolic tangent with LeCun's magic constants."""
    X *= 2/3.
    np.tanh(X, X)
    X *= 1.7159
    return X


class MLPClassifier(BaseEstimator, ClassifierMixin):
    """Multi-layer perceptron (feedforward neural network) classifier.

    Trained with gradient descent under log loss (aka the cross-entropy error
    function).

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
    momentum : float, optional
        Parameter for the momentum method. Set this somewhere between .5 and 1.
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
    def __init__(self, n_hidden, activation="tanh", alpha=0, batch_size=100,
                 learning_rate=.01, max_iter=100, momentum=.9,
                 random_state=None, shuffle=True, tol=1e-5, verbose=False):
        self.activation = activation
        self.alpha = alpha
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.max_iter = max_iter
        self.momentum = momentum
        self.n_hidden = n_hidden
        self.random_state = check_random_state(random_state)
        self.shuffle = shuffle
        self.tol = tol
        self.verbose = verbose

    def _init_fit(self, n_features, n_targets):
        n_hidden = self.n_hidden
        rng = self.random_state

        # XXX Hinton's advice is to initialize weights to be proportional to
        # sqrt(fan-in). That doesn't work on all problems, though, and it
        # blows up in the face of L2 penalty.
        #sqrt_n_features = np.sqrt(n_features)
        #sqrt_n_hidden = np.sqrt(n_hidden)

        self.coef_hidden_ = rng.uniform(-1, 1, (n_hidden, n_features))
        self.coef_output_ = rng.uniform(-1, 1, (n_targets, n_hidden))

        self.intercept_hidden_ = rng.uniform(-1, 1, n_hidden)
        self.intercept_output_ = rng.uniform(-1, 1, n_targets)

    def fit(self, X, y):
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        _, n_features = X.shape

        self._lbin = LabelBinarizer()
        Y = self._lbin.fit_transform(y)

        self._init_fit(n_features, Y.shape[1])
        backprop_sgd(self, X, Y)

        return self

    def partial_fit(self, X, y, classes):
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

        backprop_sgd(self, X, Y)

        return self

    def decision_function(self, X):
        X = atleast2d_or_csr(X)
        z_hidden = (safe_sparse_dot(X, self.coef_hidden_.T) +
                    self.intercept_hidden_)
        y_hidden = logistic(z_hidden) if self.activation == "logistic" \
                                      else _tanh(z_hidden)
        y_output = (np.dot(y_hidden, self.coef_output_.T) +
                    self.intercept_output_)
        if y_output.shape[1] == 1:
            y_output = y_output.ravel()
        return y_output

    def predict(self, X):
        scores = self.decision_function(X)
        if len(scores.shape) == 1:
            indices = (scores > 0).astype(np.int)
        else:
            indices = scores.argmax(axis=1)
        return self.classes_[indices]

    def predict_log_proba(self, X):
        scores = self.decision_function(X)
        if len(scores.shape) == 1:
            return np.log(logistic(scores))
        else:
            return log_softmax(scores)

    def predict_proba(self, X):
        scores = self.decision_function(X)
        if len(scores.shape) == 1:
            return logistic(scores)
        else:
            return softmax(scores)

    @property
    def classes_(self):
        return self._lbin.classes_
