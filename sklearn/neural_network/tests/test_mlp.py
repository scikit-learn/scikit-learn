"""
Testing for Multi-layer Perceptron module (sklearn.neural_network)
"""

# Author: Issam H. Laradji
# Licence: BSD 3 clause

import sys
import warnings

import numpy as np

from numpy.testing import assert_almost_equal, assert_array_equal

from sklearn.datasets import load_digits, load_boston
from sklearn.datasets import make_regression, make_multilabel_classification
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.metrics import roc_auc_score
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy.sparse import csr_matrix
from sklearn.utils.testing import (assert_raises, assert_greater, assert_equal,
                                   assert_false)


np.seterr(all='warn')

ACTIVATION_TYPES = ["logistic", "tanh", "relu"]

digits_dataset_multi = load_digits(n_class=3)

X_digits_multi = MinMaxScaler().fit_transform(digits_dataset_multi.data[:200])
y_digits_multi = digits_dataset_multi.target[:200]

digits_dataset_binary = load_digits(n_class=2)

X_digits_binary = MinMaxScaler().fit_transform(
    digits_dataset_binary.data[:200])
y_digits_binary = digits_dataset_binary.target[:200]

classification_datasets = [(X_digits_multi, y_digits_multi),
                           (X_digits_binary, y_digits_binary)]

boston = load_boston()

Xboston = StandardScaler().fit_transform(boston.data)[: 200]
yboston = boston.target[:200]


def test_alpha():
    # Test that larger alpha yields weights closer to zero"""
    X = X_digits_binary[:100]
    y = y_digits_binary[:100]

    alpha_vectors = []
    alpha_values = np.arange(2)
    absolute_sum = lambda x: np.sum(np.abs(x))

    for alpha in alpha_values:
        mlp = MLPClassifier(hidden_layer_sizes=10, alpha=alpha, random_state=1)
        mlp.fit(X, y)
        alpha_vectors.append(np.array([absolute_sum(mlp.coefs_[0]),
                                       absolute_sum(mlp.coefs_[1])]))

    for i in range(len(alpha_values) - 1):
        assert (alpha_vectors[i] > alpha_vectors[i + 1]).all()


def test_fit():
    # Test that the algorithm solution is equal to a worked out example."""
    X = np.array([[0.6, 0.8, 0.7]])
    y = np.array([0])
    mlp = MLPClassifier(algorithm='sgd', learning_rate_init=0.1, alpha=0.1,
                        activation='logistic', random_state=1, max_iter=1,
                        hidden_layer_sizes=2, momentum=0)
    # set weights
    mlp.coefs_ = [0] * 2
    mlp.intercepts_ = [0] * 2
    mlp.classes_ = [0, 1]
    mlp.n_outputs_ = 1
    mlp.coefs_[0] = np.array([[0.1, 0.2], [0.3, 0.1], [0.5, 0]])
    mlp.coefs_[1] = np.array([[0.1], [0.2]])
    mlp.intercepts_[0] = np.array([0.1, 0.1])
    mlp.intercepts_[1] = np.array([1.0])
    mlp._coef_grads = [] * 2
    mlp._intercept_grads = [] * 2

    mlp.label_binarizer_.y_type_ = 'binary'

    # Initialize parameters
    mlp.n_iter_ = 0
    mlp.learning_rate_ = 0.1

    # Compute the number of layers
    mlp.n_layers_ = 3

    # Pre-allocate gradient matrices
    mlp._coef_grads = [0] * (mlp.n_layers_ - 1)
    mlp._intercept_grads = [0] * (mlp.n_layers_ - 1)

    mlp.out_activation_ = 'logistic'
    mlp.t_ = 0
    mlp.best_loss_ = np.inf
    mlp.loss_curve_ = []
    mlp._no_improvement_count = 0
    mlp._intercept_velocity = [np.zeros_like(intercepts) for
                               intercepts in
                               mlp.intercepts_]
    mlp._coef_velocity = [np.zeros_like(coefs) for coefs in
                          mlp.coefs_]

    mlp.partial_fit(X, y, classes=[0, 1])
    # Manually worked out example
    # h1 = g(X1 * W_i1 + b11) = g(0.6 * 0.1 + 0.8 * 0.3 + 0.7 * 0.5 + 0.1)
    #       =  0.679178699175393
    # h2 = g(X2 * W_i2 + b12) = g(0.6 * 0.2 + 0.8 * 0.1 + 0.7 * 0 + 0.1)
    #         = 0.574442516811659
    # o1 = g(h * W2 + b21) = g(0.679 * 0.1 + 0.574 * 0.2 + 1)
    #       = 0.7654329236196236
    # d21 = -(0 - 0.765) = 0.765
    # d11 = (1 - 0.679) * 0.679 * 0.765 * 0.1 = 0.01667
    # d12 = (1 - 0.574) * 0.574 * 0.765 * 0.2 = 0.0374
    # W1grad11 = X1 * d11 + alpha * W11 = 0.6 * 0.01667 + 0.1 * 0.1 = 0.0200
    # W1grad11 = X1 * d12 + alpha * W12 = 0.6 * 0.0374 + 0.1 * 0.2 = 0.04244
    # W1grad21 = X2 * d11 + alpha * W13 = 0.8 * 0.01667 + 0.1 * 0.3 = 0.043336
    # W1grad22 = X2 * d12 + alpha * W14 = 0.8 * 0.0374 + 0.1 * 0.1 = 0.03992
    # W1grad31 = X3 * d11 + alpha * W15 = 0.6 * 0.01667 + 0.1 * 0.5 = 0.060002
    # W1grad32 = X3 * d12 + alpha * W16 = 0.6 * 0.0374 + 0.1 * 0 = 0.02244
    # W2grad1 = h1 * d21 + alpha * W21 = 0.679 * 0.765 + 0.1 * 0.1 = 0.5294
    # W2grad2 = h2 * d21 + alpha * W22 = 0.574 * 0.765 + 0.1 * 0.2 = 0.45911
    # b1grad1 = d11 = 0.01667
    # b1grad2 = d12 = 0.0374
    # b2grad = d21 = 0.765
    # W1 = W1 - eta * [W1grad11, .., W1grad32] = [[0.1, 0.2], [0.3, 0.1],
    #          [0.5, 0]] - 0.1 * [[0.0200, 0.04244], [0.043336, 0.03992],
    #          [0.060002, 0.02244]] = [[0.098, 0.195756], [0.2956664,
    #          0.096008], [0.4939998, -0.002244]]
    # W2 = W2 - eta * [W2grad1, W2grad2] = [[0.1], [0.2]] - 0.1 *
    #        [[0.5294], [0.45911]] = [[0.04706], [0.154089]]
    # b1 = b1 - eta * [b1grad1, b1grad2] = 0.1 - 0.1 * [0.01667, 0.0374]
    #         = [0.098333, 0.09626]
    # b2 = b2 - eta * b2grad = 1.0 - 0.1 * 0.765 = 0.9235
    assert_almost_equal(mlp.coefs_[0], np.array([[0.098, 0.195756],
                                                 [0.2956664, 0.096008],
                                                 [0.4939998, -0.002244]]),
                        decimal=3)
    assert_almost_equal(mlp.coefs_[1], np.array([[0.04706], [0.154089]]),
                        decimal=3)
    assert_almost_equal(mlp.intercepts_[0],
                        np.array([0.098333, 0.09626]), decimal=3)
    assert_almost_equal(mlp.intercepts_[1], np.array(0.9235), decimal=3)
    # Testing output
    #  h1 = g(X1 * W_i1 + b11) = g(0.6 * 0.098 + 0.8 * 0.2956664 +
    #               0.7 * 0.4939998 + 0.098333) = 0.677
    #  h2 = g(X2 * W_i2 + b12) = g(0.6 * 0.195756 + 0.8 * 0.096008 +
    #            0.7 * -0.002244 + 0.09626) = 0.572
    #  o1 = h * W2 + b21 = 0.677 * 0.04706 +
    #             0.572 * 0.154089 + 0.9235 = 1.043
    assert_almost_equal(mlp.decision_function(X), 1.043, decimal=3)


def test_gradient():
    # Test gradient.

    # This makes sure that the activation functions and their derivatives
    # are correct. The numerical and analytical computation of the gradient
    # should be close.
    for n_labels in [2, 3]:
        n_samples = 5
        n_features = 10
        X = np.random.random((n_samples, n_features))
        y = 1 + np.mod(np.arange(n_samples) + 1, n_labels)
        Y = LabelBinarizer().fit_transform(y)

        for activation in ACTIVATION_TYPES:
            mlp = MLPClassifier(activation=activation, hidden_layer_sizes=10,
                                algorithm='l-bfgs', alpha=1e-5,
                                learning_rate_init=0.2, max_iter=1,
                                random_state=1)
            mlp.fit(X, y)

            theta = np.hstack([l.ravel() for l in mlp.coefs_ +
                               mlp.intercepts_])

            layer_units = ([X.shape[1]] + [mlp.hidden_layer_sizes] +
                           [mlp.n_outputs_])

            activations = []
            deltas = []
            coef_grads = []
            intercept_grads = []

            activations.append(X)
            for i in range(mlp.n_layers_ - 1):
                activations.append(np.empty((X.shape[0],
                                             layer_units[i + 1])))
                deltas.append(np.empty((X.shape[0],
                                        layer_units[i + 1])))

                fan_in = layer_units[i]
                fan_out = layer_units[i + 1]
                coef_grads.append(np.empty((fan_in, fan_out)))
                intercept_grads.append(np.empty(fan_out))

            # analytically compute the gradients
            def loss_grad_fun(t):
                return mlp._loss_grad_lbfgs(t, X, Y, activations, deltas,
                                            coef_grads, intercept_grads)

            [value, grad] = loss_grad_fun(theta)
            numgrad = np.zeros(np.size(theta))
            n = np.size(theta, 0)
            E = np.eye(n)
            epsilon = 1e-5
            # numerically compute the gradients
            for i in range(n):
                dtheta = E[:, i] * epsilon
                numgrad[i] = ((loss_grad_fun(theta + dtheta)[0] -
                              loss_grad_fun(theta - dtheta)[0]) /
                              (epsilon * 2.0))
            assert_almost_equal(numgrad, grad)


def test_lbfgs_classification():
    # Test lbfgs on classification.
    # It should achieve a score higher than 0.95 for the binary and multi-class
    # versions of the digits dataset.
    for X, y in classification_datasets:
        X_train = X[:150]
        y_train = y[:150]
        X_test = X[150:]

        expected_shape_dtype = (X_test.shape[0], y_train.dtype.kind)

        for activation in ACTIVATION_TYPES:
            mlp = MLPClassifier(algorithm='l-bfgs', hidden_layer_sizes=50,
                                max_iter=150, shuffle=True, random_state=1,
                                activation=activation)
            mlp.fit(X_train, y_train)
            y_predict = mlp.predict(X_test)
            assert_greater(mlp.score(X_train, y_train), 0.95)
            assert_equal((y_predict.shape[0], y_predict.dtype.kind),
                         expected_shape_dtype)


def test_lbfgs_regression():
    # Test lbfgs on the boston dataset, a regression problems."""
    X = Xboston
    y = yboston
    for activation in ACTIVATION_TYPES:
        mlp = MLPRegressor(algorithm='l-bfgs', hidden_layer_sizes=50,
                           max_iter=150, shuffle=True, random_state=1,
                           activation=activation)
        mlp.fit(X, y)
        assert_greater(mlp.score(X, y), 0.95)


def test_learning_rate_warmstart():
    # Tests that warm_start reuses past solution."""
    X = [[3, 2], [1, 6], [5, 6], [-2, -4]]
    y = [1, 1, 1, 0]
    for learning_rate in ["invscaling", "constant"]:
        mlp = MLPClassifier(algorithm='sgd', hidden_layer_sizes=4,
                            learning_rate=learning_rate, max_iter=1,
                            power_t=0.25, warm_start=True)
        mlp.fit(X, y)
        prev_eta = mlp._optimizer.learning_rate
        mlp.fit(X, y)
        post_eta = mlp._optimizer.learning_rate

        if learning_rate == 'constant':
            assert_equal(prev_eta, post_eta)
        elif learning_rate == 'invscaling':
            assert_equal(mlp.learning_rate_init / pow(8 + 1, mlp.power_t),
                         post_eta)


def test_multilabel_classification():
    # Test that multi-label classification works as expected."""
    # test fit method
    X, y = make_multilabel_classification(n_samples=50, random_state=0,
                                          return_indicator=True)
    mlp = MLPClassifier(algorithm='l-bfgs', hidden_layer_sizes=50, alpha=1e-5,
                        max_iter=150, random_state=0, activation='logistic',
                        learning_rate_init=0.2)
    mlp.fit(X, y)
    assert_equal(mlp.score(X, y), 1)

    # test partial fit method
    mlp = MLPClassifier(algorithm='sgd', hidden_layer_sizes=50, max_iter=150,
                        random_state=0, activation='logistic', alpha=1e-5,
                        learning_rate_init=0.2)
    for i in range(100):
        mlp.partial_fit(X, y, classes=[0, 1, 2, 3, 4])
    assert_greater(mlp.score(X, y), 0.9)


def test_multioutput_regression():
    # Test that multi-output regression works as expected"""
    X, y = make_regression(n_samples=200, n_targets=5)
    mlp = MLPRegressor(algorithm='l-bfgs', hidden_layer_sizes=50, max_iter=200,
                       random_state=1)
    mlp.fit(X, y)
    assert_greater(mlp.score(X, y), 0.9)


def test_partial_fit_classes_error():
    # Tests that passing different classes to partial_fit raises an error"""
    X = [[3, 2]]
    y = [0]
    clf = MLPClassifier(algorithm='sgd')
    clf.partial_fit(X, y, classes=[0, 1])
    assert_raises(ValueError, clf.partial_fit, X, y, classes=[1, 2])


def test_partial_fit_classification():
    # Test partial_fit on classification.
    # `partial_fit` should yield the same results as 'fit'for binary and
    # multi-class classification.
    for X, y in classification_datasets:
        X = X
        y = y
        mlp = MLPClassifier(algorithm='sgd', max_iter=100, random_state=1,
                            tol=0, alpha=1e-5, learning_rate_init=0.2)

        mlp.fit(X, y)
        pred1 = mlp.predict(X)
        mlp = MLPClassifier(algorithm='sgd', random_state=1, alpha=1e-5,
                            learning_rate_init=0.2)
        for i in range(100):
            mlp.partial_fit(X, y, classes=np.unique(y))
        pred2 = mlp.predict(X)
        assert_array_equal(pred1, pred2)
        assert_greater(mlp.score(X, y), 0.95)


def test_partial_fit_regression():
    # Test partial_fit on regression.
    # `partial_fit` should yield the same results as 'fit' for regression.
    X = Xboston
    y = yboston

    for momentum in [0, .9]:
        mlp = MLPRegressor(algorithm='sgd', max_iter=100, activation='relu',
                           random_state=1, learning_rate_init=0.01,
                           batch_size=X.shape[0], momentum=momentum)
        with warnings.catch_warnings(record=True):
            # catch convergence warning
            mlp.fit(X, y)
        pred1 = mlp.predict(X)
        mlp = MLPRegressor(algorithm='sgd', activation='relu',
                           learning_rate_init=0.01, random_state=1,
                           batch_size=X.shape[0], momentum=momentum)
        for i in range(100):
            mlp.partial_fit(X, y)

        pred2 = mlp.predict(X)
        assert_almost_equal(pred1, pred2, decimal=2)
        score = mlp.score(X, y)
        assert_greater(score, 0.75)


def test_partial_fit_errors():
    # Test partial_fit error handling."""
    X = [[3, 2], [1, 6]]
    y = [1, 0]

    # no classes passed
    assert_raises(ValueError,
                  MLPClassifier(
                      algorithm='sgd').partial_fit,
                  X, y,
                  classes=[2])

    # l-bfgs doesn't support partial_fit
    assert_false(hasattr(MLPClassifier(algorithm='l-bfgs'), 'partial_fit'))


def test_params_errors():
    # Test that invalid parameters raise value error"""
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MLPClassifier

    assert_raises(ValueError, clf(hidden_layer_sizes=-1).fit, X, y)
    assert_raises(ValueError, clf(max_iter=-1).fit, X, y)
    assert_raises(ValueError, clf(shuffle='true').fit, X, y)
    assert_raises(ValueError, clf(alpha=-1).fit, X, y)
    assert_raises(ValueError, clf(learning_rate_init=-1).fit, X, y)

    assert_raises(ValueError, clf(algorithm='hadoken').fit, X, y)
    assert_raises(ValueError, clf(learning_rate='converge').fit, X, y)
    assert_raises(ValueError, clf(activation='cloak').fit, X, y)


def test_predict_proba_binary():
    # Test that predict_proba works as expected for binary class."""
    X = X_digits_binary[:50]
    y = y_digits_binary[:50]

    clf = MLPClassifier(hidden_layer_sizes=5)
    clf.fit(X, y)
    y_proba = clf.predict_proba(X)
    y_log_proba = clf.predict_log_proba(X)

    (n_samples, n_classes) = y.shape[0], 2

    proba_max = y_proba.argmax(axis=1)
    proba_log_max = y_log_proba.argmax(axis=1)

    assert_equal(y_proba.shape, (n_samples, n_classes))
    assert_array_equal(proba_max, proba_log_max)
    assert_array_equal(y_log_proba, np.log(y_proba))

    assert_equal(roc_auc_score(y, y_proba[:, 1]), 1.0)


def test_predict_proba_multi():
    # Test that predict_proba works as expected for multi class."""
    X = X_digits_multi[:10]
    y = y_digits_multi[:10]

    clf = MLPClassifier(hidden_layer_sizes=5)
    clf.fit(X, y)
    y_proba = clf.predict_proba(X)
    y_log_proba = clf.predict_log_proba(X)

    (n_samples, n_classes) = y.shape[0], np.unique(y).size

    proba_max = y_proba.argmax(axis=1)
    proba_log_max = y_log_proba.argmax(axis=1)

    assert_equal(y_proba.shape, (n_samples, n_classes))
    assert_array_equal(proba_max, proba_log_max)
    assert_array_equal(y_log_proba, np.log(y_proba))


def test_sparse_matrices():
    # Test that sparse and dense input matrices output the same results."""
    X = X_digits_binary[:50]
    y = y_digits_binary[:50]
    X_sparse = csr_matrix(X)
    mlp = MLPClassifier(random_state=1, hidden_layer_sizes=15)
    mlp.fit(X, y)
    pred1 = mlp.decision_function(X)
    mlp.fit(X_sparse, y)
    pred2 = mlp.decision_function(X_sparse)
    assert_almost_equal(pred1, pred2)
    pred1 = mlp.predict(X)
    pred2 = mlp.predict(X_sparse)
    assert_array_equal(pred1, pred2)


def test_tolerance():
    # Test tolerance.
    # It should force the algorithm to exit the loop when it converges.
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MLPClassifier(tol=0.5, max_iter=3000, algorithm='sgd', verbose=10)
    clf.fit(X, y)
    assert_greater(clf.max_iter, clf.n_iter_)


def test_verbose_sgd():
    # Test verbose.
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MLPClassifier(algorithm='sgd', max_iter=2, verbose=10,
                        hidden_layer_sizes=2)
    old_stdout = sys.stdout
    sys.stdout = output = StringIO()

    clf.fit(X, y)
    clf.partial_fit(X, y)

    sys.stdout = old_stdout
    assert 'Iteration' in output.getvalue()


def test_early_stopping():
    X = X_digits_binary[:100]
    y = y_digits_binary[:100]
    tol = 0.2
    clf = MLPClassifier(tol=tol, max_iter=3000, algorithm='sgd',
                        early_stopping=True)
    clf.fit(X, y)
    assert_greater(clf.max_iter, clf.n_iter_)

    valid_scores = clf.validation_scores_
    best_valid_score = clf.best_validation_score_
    assert_equal(max(valid_scores), best_valid_score)
    assert_greater(best_valid_score + tol, valid_scores[-2])
    assert_greater(best_valid_score + tol, valid_scores[-1])


def test_adaptive_learning_rate():
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MLPClassifier(tol=0.5, max_iter=3000, algorithm='sgd',
                        learning_rate='adaptive', verbose=10)
    clf.fit(X, y)
    assert_greater(clf.max_iter, clf.n_iter_)
    assert_greater(1e-6, clf._optimizer.learning_rate)
