"""
Testing for Multi-layer Perceptron module (sklearn.neural_network)
"""

# Author: Issam H. Laradji
# Licence: BSD 3 clause

import sys

import numpy as np

from numpy.testing import assert_almost_equal, assert_array_equal

from sklearn.datasets import load_digits, load_boston
from sklearn.datasets import make_regression, make_multilabel_classification
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.metrics import auc_score
from sklearn.neural_network import MultilayerPerceptronClassifier
from sklearn.neural_network import MultilayerPerceptronRegressor
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import StandardScaler
from scipy.sparse import csr_matrix
from sklearn.utils.testing import assert_raises, assert_greater, assert_equal


np.seterr(all='warn')

LEARNING_RATE_TYPES = ["constant", "invscaling"]

ACTIVATION_TYPES = ["logistic", "tanh"]

digits_dataset_multi = load_digits(n_class=3)

Xdigits_multi = digits_dataset_multi.data[:200]
ydigits_multi = digits_dataset_multi.target[:200]
Xdigits_multi -= Xdigits_multi.min()
Xdigits_multi /= Xdigits_multi.max()

digits_dataset_binary = load_digits(n_class=2)

Xdigits_binary = digits_dataset_binary.data[:200]
ydigits_binary = digits_dataset_binary.target[:200]
Xdigits_binary -= Xdigits_binary.min()
Xdigits_binary /= Xdigits_binary.max()

classification_datasets = [(Xdigits_multi, ydigits_multi),
                           (Xdigits_binary, ydigits_binary)]

boston = load_boston()

Xboston = StandardScaler().fit_transform(boston.data)[: 200]
yboston = boston.target[:200]


def test_alpha():
    """Tests that larger alpha  yields weights closer  to zero"""
    X = Xdigits_binary[:100]
    y = ydigits_binary[:100]

    alpha_vectors = []
    alpha_values = np.arange(2)
    absolute_sum = lambda x: np.sum(np.abs(x))

    for alpha in alpha_values:
        mlp = MultilayerPerceptronClassifier(
            n_hidden=10,
            alpha=alpha, random_state=1)
        mlp.fit(X, y)
        alpha_vectors.append(np.array([absolute_sum(mlp.coef_hidden_),
                                       absolute_sum(mlp.coef_output_)]))

    for i in range(len(alpha_values) - 1):
        assert (alpha_vectors[i] > alpha_vectors[i + 1]).all()


def test_fit():
    """
    Tests that algorithm solution is equal to
    a worked out example
    """
    X = np.array([[0.6, 0.8, 0.7]])
    y = np.array([0])
    mlp = MultilayerPerceptronClassifier(algorithm='sgd',
                                         eta0=0.1, alpha=0.1,
                                         activation='logistic', random_state=1,
                                         max_iter=1, n_hidden=2)
    # set weights
    mlp.classes_ = [0, 1]
    mlp.coef_hidden_ = np.array([[0.1, 0.2], [0.3, 0.1], [0.5, 0]])
    mlp.coef_output_ = np.array([[0.1], [0.2]])
    mlp.intercept_hidden_ = np.array([0.1, 0.1])
    mlp.intercept_output_ = np.array([1.0])
    # initialize
    mlp._init_param()
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
    assert_almost_equal(mlp.coef_hidden_, np.array(
        [[0.098, 0.195756], [0.2956664, 0.096008],
        [0.4939998, -0.002244]]), decimal=3)
    assert_almost_equal(
        mlp.coef_output_,
        np.array([[0.04706],
                  [0.154089]]),
        decimal=3)
    assert_almost_equal(
        mlp.intercept_hidden_,
        np.array([0.098333,
                  0.09626]),
        decimal=3)
    assert_almost_equal(mlp.intercept_output_, np.array(0.9235), decimal=3)
    # Testing output
    #  h1 = g(X1 * W_i1 + b11) = g(0.6 * 0.098 + 0.8 * 0.2956664 +
    #               0.7 * 0.4939998 + 0.098333) = 0.677
    #  h2 = g(X2 * W_i2 + b12) = g(0.6 * 0.195756 + 0.8 * 0.096008 +
    #            0.7 * -0.002244 + 0.09626) = 0.572
    #  o1 = g(h * W2 + b21) = g(0.677 * 0.04706 +
    #             0.572 * 0.154089 + 0.9235) = 1.043
    assert_almost_equal(mlp.decision_function(X), 1.043, decimal=3)


def test_gradient():
    """
    Tests that numerical and analytical computation
    of the gradient are close
    """
    n_labels = 2
    n_samples = 5
    n_features = 2
    X = np.random.random((n_samples, n_features))
    y = 1 + np.mod(np.arange(n_samples) + 1, n_labels)
    Y = LabelBinarizer().fit_transform(y)
    n_outputs = Y.shape[1]
    mlp = MultilayerPerceptronClassifier(alpha=0.5, n_hidden=2)
    mlp.classes_ = np.unique(Y)
    mlp.n_features = n_features
    mlp.n_outputs = n_outputs
    mlp._init_fit()
    mlp._init_param()
    theta = mlp._pack(mlp.coef_hidden_, mlp.coef_output_,
                      mlp.intercept_hidden_, mlp.intercept_output_)
    a_hidden, a_output, delta_o = mlp._preallocate_memory(n_samples)
    # analytically compute the gradients
    J = lambda t: mlp._cost_grad(
        t,
        X,
        Y,
        n_samples,
        a_hidden,
        a_output,
        delta_o)
    [value, grad] = J(theta)
    numgrad = np.zeros(np.size(theta))
    n = np.size(theta, 0)
    E = np.eye(n)
    epsilon = 1e-4
    # numerically compute the gradients
    for i in range(n):
        dtheta = E[:, i] * epsilon
        numgrad[i] = (
            J(theta + dtheta)[0] - J(theta - dtheta)[0]) / (epsilon * 2.0)
    assert_almost_equal(numgrad, grad)


def test_lbfgs_classification():
    """
    Tests that L-bfgs achieves score higher than 0.95 for binary-
    and multi-classification digits datasets
    """
    for X, y in classification_datasets:
        X_train = X[:150]
        y_train = y[:150]
        X_test = X[150:]

        expected_shape_dtype = (X_test.shape[0], y_train.dtype.kind)

        for activation in ACTIVATION_TYPES:
            mlp = MultilayerPerceptronClassifier(
                algorithm='l-bfgs',
                n_hidden=50,
                max_iter=150,
                shuffle=True,
                random_state=1,
                activation=activation)
            mlp.fit(X_train, y_train)
            y_predict = mlp.predict(X_test)

            assert_greater(mlp.score(X_train, y_train), 0.9)
            assert_equal(
                (y_predict.shape[0],
                 y_predict.dtype.kind),
                expected_shape_dtype)


def test_lbfgs_regression():
    """
    Tests that L-bfgs achieves score higher than 0.95 for the
    boston dataset
    """
    X = Xboston
    y = yboston
    for activation in ACTIVATION_TYPES:
        mlp = MultilayerPerceptronRegressor(
            algorithm='l-bfgs',
            n_hidden=50,
            max_iter=150,
            shuffle=True,
            random_state=1,
            activation=activation)
        mlp.fit(X, y)
        assert_greater(mlp.score(X, y), 0.95)


def test_learning_rate_warmstart():
    """
    Tests that warm_start reuses past solution.

    Learning rates should behave as expected.
    """
    X = [[3, 2], [1, 6], [5, 6], [-2, -4]]
    y = [1, 1, 1, 0]
    for learning_rate in LEARNING_RATE_TYPES:
        mlp = MultilayerPerceptronClassifier(
            algorithm='sgd', n_hidden=4,
            learning_rate=learning_rate,
            max_iter=1, power_t=0.25,
            warm_start=True)
        mlp.fit(X, y)
        prev_eta = mlp.eta_
        mlp.fit(X, y)
        post_eta = mlp.eta_
        if learning_rate == 'constant':
            assert prev_eta == post_eta
        elif learning_rate == 'invscaling':
            assert post_eta == prev_eta / pow(mlp.t_ - 1, mlp.power_t)


def test_multilabel_classification():
    """
    Tests that multi-label classification works as expected
    """
    # test fit method
    X, y = make_multilabel_classification(n_samples=50, random_state=0)
    mlp = MultilayerPerceptronClassifier(
        algorithm='l-bfgs',
        n_hidden=50,
        max_iter=150,
        random_state=1,
        activation='logistic')
    mlp.fit(X, y)
    assert_equal(mlp.score(X, y), 1)

    # test partial fit method
    mlp = MultilayerPerceptronClassifier(
        algorithm='sgd',
        n_hidden=50,
        max_iter=150,
        random_state=1,
        activation='logistic')
    y = np.atleast_1d(y)
    classes_ = np.unique(y)
    for i in xrange(100):
            mlp.partial_fit(X, y, classes=classes_)
    assert_greater(mlp.score(X, y), 0.95)


def test_multioutput_regression():
    """
    Tests that multi-output regression works as expected
    """
    X, y = make_regression(n_samples=200, n_targets=5)
    mlp = MultilayerPerceptronRegressor(
        algorithm='l-bfgs',
        n_hidden=50,
        max_iter=150,
        random_state=1,
        activation='logistic')
    mlp.fit(X, y)
    assert_greater(mlp.score(X, y), 0.95)


def test_partial_fit_classes_error():
    """Tests that passing different classes to partial_fit raises an error"""
    X = [3, 2]
    y = [0]
    clf = MultilayerPerceptronClassifier(algorithm='sgd')
    clf.partial_fit(X, y, classes=[0, 1])
    assert_raises(ValueError, clf.partial_fit, X, y, classes=[0, 1, 2])


def test_partial_fit_classification():
    """
    Tests that partial_fit yields same results as 'fit'
    for binary- and multi-class classification  with
    different activations functions'partial fit' should yield
    """
    for X, y in classification_datasets:
        X = X
        y = y
        mlp = MultilayerPerceptronClassifier(
            algorithm='sgd',
            max_iter=100,
            random_state=1,
            batch_size=X.shape[0])
        mlp.fit(X, y)
        pred1 = mlp.predict(X)
        mlp = MultilayerPerceptronClassifier(algorithm='sgd', random_state=1)
        for i in xrange(100):
            mlp.partial_fit(X, y, classes=np.unique(y))
        pred2 = mlp.predict(X)
        assert_array_equal(pred1, pred2)
        assert_greater(mlp.score(X, y), 0.95)


def test_partial_fit_regression():
    """
    Tests that partial_fit yields same results as 'fit'
    for regression  with different activations functions
    """
    X = Xboston
    y = yboston
    for activation in ACTIVATION_TYPES:
        mlp = MultilayerPerceptronRegressor(
            algorithm='sgd',
            max_iter=100,
            random_state=1,
            activation=activation,
            batch_size=X.shape[0])
        mlp.fit(X, y)
        pred1 = mlp.predict(X)
        mlp = MultilayerPerceptronRegressor(
            algorithm='sgd',
            activation=activation,
            random_state=1)
        for i in xrange(100):
            mlp.partial_fit(X, y)
        pred2 = mlp.predict(X)
        assert_almost_equal(pred1, pred2, decimal=2)
        print mlp.score(X, y)
        assert_greater(mlp.score(X, y), 0.9)


def test_partial_fit_errors():
    """
    Tests that partial_fit raises value errors when given invalid parameters.
    """
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MultilayerPerceptronClassifier

    # no classes passed
    assert_raises(ValueError, clf(algorithm='sgd').partial_fit, X, y)
    assert_raises(
        ValueError,
        clf(algorithm='l-bfgs').partial_fit,
        X, y,
        classes=[0, 1])


def test_params_errors():
    """Tests that invalid parameters raise value error"""
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MultilayerPerceptronClassifier

    assert_raises(ValueError, clf(n_hidden=-1).fit, X, y)
    assert_raises(ValueError, clf(max_iter=-1).fit, X, y)
    assert_raises(ValueError, clf(shuffle='true').fit, X, y)
    assert_raises(ValueError, clf(alpha=-1).fit, X, y)
    assert_raises(ValueError, clf(eta0=-1).fit, X, y)

    assert_raises(ValueError, clf(algorithm='hadoken').fit, X, y)
    assert_raises(ValueError, clf(learning_rate='converge').fit, X, y)
    assert_raises(ValueError, clf(activation='cloak').fit, X, y)


def test_predict_proba_binary():
    """
    Tests that predict_proba works as expected for binary class.
    """
    X = Xdigits_binary[:50]
    y = ydigits_binary[:50]

    clf = MultilayerPerceptronClassifier(n_hidden=5)
    clf.fit(X, y)
    y_proba = clf.predict_proba(X)
    y_log_proba = clf.predict_log_proba(X)

    (n_samples, n_classes) = y.shape[0], 2

    proba_max = y_proba.argmax(axis=1)
    proba_log_max = y_log_proba.argmax(axis=1)

    assert_equal(y_proba.shape, (n_samples, n_classes))
    assert_array_equal(proba_max, proba_log_max)
    assert_array_equal(y_log_proba, np.log(y_proba))

    assert_equal(auc_score(y, y_proba[:, 1]), 1.0)


def test_predict_proba_multi():
    """
    Tests that predict_proba works as expected for multi class.
    """
    X = Xdigits_multi[:10]
    y = ydigits_multi[:10]

    clf = MultilayerPerceptronClassifier(n_hidden=5)
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
    """
    Tests that sparse and dense input matrices
    yield equal output
    """
    X = Xdigits_binary[:50]
    y = ydigits_binary[:50]
    X_sparse = csr_matrix(X)
    mlp = MultilayerPerceptronClassifier(random_state=1, n_hidden=15)
    mlp.fit(X, y)
    pred1 = mlp.decision_function(X)
    mlp.fit(X_sparse, y)
    pred2 = mlp.decision_function(X_sparse)
    assert_almost_equal(pred1, pred2)
    pred1 = mlp.predict(X)
    pred2 = mlp.predict(X_sparse)
    assert_array_equal(pred1, pred2)


def test_tolerance():
    """
    Tests that tolerance exits the algorithm before it
    iterates for max_iter times.
    """
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MultilayerPerceptronClassifier(
        tol=0.5,
        max_iter=3000,
        algorithm='sgd')
    clf.fit(X, y)
    assert clf.t_ != clf.max_iter


def test_verbose_sgd():
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = MultilayerPerceptronClassifier(
        algorithm='sgd',
        max_iter=2,
        verbose=10,
        n_hidden=2)
    old_stdout = sys.stdout
    sys.stdout = output = StringIO()

    clf.fit(X, y)
    clf.partial_fit(X, y)

    sys.stdout = old_stdout
    assert 'Iteration' in output.getvalue()
