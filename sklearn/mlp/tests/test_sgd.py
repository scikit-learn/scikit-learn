import numpy as np

from ..mlp_fast import forward, backward, LogSig, Tanh, CrossEntropyLoss, SquaredLoss
from ..classes import MLPClassifier
from ... import datasets, preprocessing

from nose.tools import assert_true
from scipy.optimize import check_grad


def test_cross_entropy():
    X, y = datasets.make_classification()

    X = preprocessing.scale(X)

    classifier = MLPClassifier(n_hidden=10, lr=0.3, lr_moment=0.3,
        batch_size=100, loss_function='cross-entropy')
    classifier.fit(X, y,
        max_epochs=500)

    assert_true(np.mean(classifier.predict(X) == y) > .95)


def test_numbers():
    digits = datasets.load_digits()
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))

    data = preprocessing.scale(data)

    classifier = MLPClassifier(n_hidden=10, lr=0.3, lr_moment=0.3, batch_size=100)
    classifier.fit(data, digits.target,
        max_epochs=500)

    assert_true(np.mean(classifier.predict(data) == digits.target) > .95)


def test_iris():
    iris = datasets.load_iris()
    data, target = iris.data, iris.target
    data = preprocessing.scale(data)

    classifier = MLPClassifier(n_hidden=10, lr=0.3, lr_moment=0.3, batch_size=1)
    classifier.fit(data, target, max_epochs=300)

    assert_true(np.mean(classifier.predict(data) == target) > .95)


def test_gradient():
    n_features = 2
    n_hidden = 2
    n_outs = 1
    batch_size = 1

    X = np.random.normal(size=(batch_size, n_features))
    y = np.random.normal(size=(batch_size, n_outs))
    x_hidden = np.empty((batch_size, n_hidden))
    x_output = np.empty((batch_size, n_outs))

    weights_hidden = np.random.normal(size=(n_features, n_hidden)) \
        / np.sqrt(n_features)
    bias_hidden = np.zeros(n_hidden)
    weights_output = np.random.normal(size=(n_hidden, n_outs)) \
        / np.sqrt(n_hidden)
    bias_output = np.zeros(n_outs)

    delta_h = np.empty((batch_size, n_hidden))
    dx_hidden = np.empty((batch_size, n_hidden))
    delta_o = np.empty((batch_size, n_outs))
    dx_output = np.empty((batch_size, n_outs))
    weights_moment_o = np.zeros((n_hidden, n_outs))
    weights_moment_h = np.zeros((n_features, n_hidden))

    def function(w0, output, loss):
        new_weights_h = w0[:n_features * n_hidden].reshape((n_features, n_hidden))
        new_weights_o = w0[n_features * n_hidden:].reshape((n_hidden, n_outs))

        forward(X, new_weights_h, bias_hidden, new_weights_o, bias_output,
            x_hidden, x_output, output, Tanh())

        out = np.empty((batch_size, n_outs))
        loss.loss(y, x_output, out)

        return out

    def gradient(w0, output, loss):
        new_weights_h = np.array(w0[:n_features * n_hidden].reshape((n_features, n_hidden)))
        new_weights_o = np.array(w0[n_features * n_hidden:].reshape((n_hidden, n_outs)))

        new_bias_o = np.array(bias_output)
        new_bias_h = np.array(bias_hidden)

        forward(X, new_weights_h, bias_hidden, new_weights_o, bias_output,
            x_hidden, x_output, output, Tanh())
        backward(X, x_output, x_hidden, y, new_weights_o, new_bias_o,
            new_weights_h, new_bias_h, delta_o, delta_h, dx_output,
            dx_hidden, weights_moment_o, weights_moment_h, loss,
            output, Tanh(), 1, 0)

        g = np.dot(X.T, -delta_h).reshape(-1)
        g = np.append(g, np.dot(x_hidden.T, -delta_o).reshape(-1))

        return g

    w = weights_hidden.reshape(-1)
    w = np.append(w, weights_output.reshape(-1))

    assert_true(check_grad(function, gradient, w, Tanh(), SquaredLoss()) < 1e-5)
    assert_true(check_grad(function, gradient, w, LogSig(), CrossEntropyLoss()) < 1e-5)
