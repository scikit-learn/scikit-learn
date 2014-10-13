"""
Testing for Extreme Learning Machines module (sklearn.neural_network)
"""

# Author: Issam H. Laradji
# Licence: BSD 3 clause

import sys

import numpy as np

from numpy.testing import assert_almost_equal, assert_array_equal
from scipy.sparse import csr_matrix
from itertools import product

from sklearn import cross_validation
from sklearn.datasets import load_digits, load_boston
from sklearn.datasets import make_regression, make_multilabel_classification
from sklearn.ensemble import AdaBoostClassifier
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.metrics import roc_auc_score
from sklearn.neural_network import ELMClassifier
from sklearn.neural_network import ELMRegressor
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.utils import gen_batches
from sklearn.utils.testing import assert_raises, assert_greater, assert_equal
from sklearn.utils.testing import assert_not_equal


np.seterr(all='warn')

random_state = 1

ACTIVATION_TYPES = ["logistic", "tanh", "relu"]
CLASSIFICATION_TYPES = ["binary", "multi-class"]

digits_dataset_multi = load_digits(n_class=3)

Xdigits_multi = MinMaxScaler().fit_transform(digits_dataset_multi.data[:200])
ydigits_multi = digits_dataset_multi.target[:200]

digits_dataset_binary = load_digits(n_class=2)

Xdigits_binary = MinMaxScaler().fit_transform(digits_dataset_binary.data[:200])
ydigits_binary = digits_dataset_binary.target[:200]


classification_datasets = {'binary': (Xdigits_binary, ydigits_binary),
                           'multi-class': (Xdigits_multi, ydigits_multi)}

boston = load_boston()

Xboston = StandardScaler().fit_transform(boston.data)[:200]
yboston = boston.target[:200]


def check_elm(elm, type_, dataset, n_samples, score):
    if type_ == 'classification':
        X, y = classification_datasets[dataset]
    elif type_ == 'regression':
        X, y = Xboston, yboston

    X_train, y_train = X[:n_samples], y[:n_samples]
    X_test = X[n_samples:]

    expected_shape = X_test.shape[0]
    expected_dtype = y_train.dtype.kind

    for activation in ACTIVATION_TYPES:
        elm.fit(X_train, y_train)

        y_predict = elm.predict(X_test)
        assert_greater(elm.score(X_train, y_train), score)

        assert_equal(y_predict.shape[0], expected_shape)
        assert_equal(y_predict.dtype.kind, expected_dtype)


def test_classification():
    """Test ELMClassifier.

    It should score higher than 0.95 for binary and multi-class
    classification digits datasets.
    """
    for name, activation in product(classification_datasets, ACTIVATION_TYPES):
        elm = ELMClassifier(n_hidden=50, activation=activation,
                            weight_scale=10, random_state=random_state)
        check_elm(elm, 'classification', name, 150, 0.95)


def test_regression():
    """Test ELMRegressor.

    It should achieve an R^2 score higher than 0.9 for the boston dataset.
    """
    for activation in ACTIVATION_TYPES:
        elm = ELMRegressor(activation=activation, C=100)
        check_elm(elm, 'regression', None, 50, 0.85)


def test_multilabel_classification():
    """Test that multi-label classification works as expected."""
    # test fit method
    X, y = make_multilabel_classification(n_samples=50, random_state=0,
                                          return_indicator=True)
    elm = ELMClassifier(weight_scale=100)
    elm.fit(X, y)
    assert_greater(elm.score(X, y), 0.95)


def test_multioutput_regression():
    """Test whether multi-output regression works as expected."""
    X, y = make_regression(n_samples=200, n_targets=5,
                           random_state=random_state)
    for activation in ACTIVATION_TYPES:
        elm = ELMRegressor(n_hidden=300, activation=activation,
                           random_state=random_state)
        elm.fit(X, y)
        assert_greater(elm.score(X, y), 0.95)


def test_non_uniform_weights_toy_edge():
    X = [[1, 0], [1, 0], [1, 0], [0, 1]]
    y = [0, 0, 1, 0]
    # ignore the first 2 training samples by setting their weight to 0
    sample_weight = [0, 0, 1, 1]
    for activation in ACTIVATION_TYPES:
        for elm in [ELMClassifier(), ELMRegressor(C=100)]:
            elm.fit(X, y, sample_weight=sample_weight)
            assert_greater(elm.predict([[1, 0]])[0], 0.5)


def test_overfitting():
    """Larger number of hidden neurons should increase training score."""
    X, y = classification_datasets['multi-class']

    for activation in ACTIVATION_TYPES:
        elm = ELMClassifier(n_hidden=5, activation=activation,
                            random_state=random_state)
        elm.fit(X, y)
        score_5_n_hidden = elm.score(X, y)

        elm = ELMClassifier(n_hidden=15, activation=activation,
                            random_state=random_state)
        elm.fit(X, y)
        score_15_n_hidden = elm.score(X, y)

        assert_greater(score_15_n_hidden, score_5_n_hidden)


def test_params_errors():
    """Test whether invalid parameters raise value error."""
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    estimator_type = ELMClassifier

    assert_raises(ValueError, estimator_type(n_hidden=-1).fit, X, y)
    assert_raises(ValueError, estimator_type(activation='ghost').fit, X, y)
    assert_raises(ValueError, estimator_type(C=-1).fit, X, y)


def test_partial_fit_classes_error():
    """Test that passing different classes to partial_fit raises an error."""
    X = [3, 2]
    y = [0]
    elm = ELMClassifier()
    # different classes passed
    assert_raises(ValueError, elm.partial_fit, X, y, classes=[1, 2])


def test_partial_fit_classification():
    """Test partial_fit for classification.

    It should output the same results as 'fit' for binary and
    multi-class classification.
    """
    for X, y in classification_datasets.values():
        batch_size = 100
        n_samples = X.shape[0]

        elm_fit = ELMClassifier(random_state=random_state,
                                batch_size=batch_size)
        elm_partial_fit = ELMClassifier(random_state=random_state)

        elm_fit.fit(X, y)
        for batch_slice in gen_batches(n_samples, batch_size):
            elm_partial_fit.partial_fit(X[batch_slice], y[batch_slice],
                                        classes=np.unique(y))

        pred1 = elm_fit.predict(X)
        pred2 = elm_partial_fit.predict(X)

        assert_array_equal(pred1, pred2)
        assert_greater(elm_fit.score(X, y), 0.95)
        assert_greater(elm_partial_fit.score(X, y), 0.95)


def test_partial_fit_regression():
    """Test partial_fit for regression.

    It should output the same results as 'fit' for regression on
    different activations functions.
    """
    X = Xboston
    y = yboston
    batch_size = 100
    n_samples = X.shape[0]

    for activation in ACTIVATION_TYPES:
        elm_fit = ELMRegressor(random_state=random_state, C=100,
                               activation=activation, batch_size=batch_size)

        elm_partial_fit = ELMRegressor(activation=activation,
                                       C=100,
                                       random_state=random_state,
                                       batch_size=batch_size)

        elm_fit.fit(X, y)
        for batch_slice in gen_batches(n_samples, batch_size):
            elm_partial_fit.partial_fit(X[batch_slice], y[batch_slice])

        pred1 = elm_fit.predict(X)
        pred2 = elm_partial_fit.predict(X)

        assert_almost_equal(pred1, pred2, decimal=2)
        assert_greater(elm_fit.score(X, y), 0.85)
        assert_greater(elm_partial_fit.score(X, y), 0.85)


def test_predict_proba_binary():
    """Test whether predict_proba works as expected for binary class."""
    X = Xdigits_binary[:50]
    y = ydigits_binary[:50]

    clf = ELMClassifier(n_hidden=10)
    clf.fit(X, y)
    y_proba = clf.predict_proba(X)

    (n_samples, n_classes) = y.shape[0], 2

    assert_equal(y_proba.shape, (n_samples, n_classes))
    assert_greater(roc_auc_score(y, y_proba[:, 1]), 0.95)


def test_predict_proba_multi():
    """Test whether predict_proba works as expected for multi class."""
    X = Xdigits_multi[:10]
    y = ydigits_multi[:10]

    clf = ELMClassifier(n_hidden=5)
    clf.fit(X, y)
    y_proba = clf.predict_proba(X)

    (n_samples, n_classes) = y.shape[0], np.unique(y).size

    assert_equal(y_proba.shape, (n_samples, n_classes))


def test_recursive_and_standard():
    """Test that recursive lsqr return the same result as standard lsqr."""
    batch_size = 50
    for dataset, class_weight in product(classification_datasets.values(),
                                         [None, 'auto']):
        X, y = dataset
        elm_standard = ELMClassifier(class_weight=class_weight,
                                     random_state=random_state)
        elm_recursive = ELMClassifier(class_weight=class_weight,
                                      random_state=random_state,
                                      batch_size=batch_size)
        elm_standard.fit(X, y)
        elm_recursive.fit(X, y)

        pred1 = elm_standard.predict(X)
        pred2 = elm_recursive.predict(X)

        assert_array_equal(pred1, pred2)
        assert_greater(elm_standard.score(X, y), 0.95)


def test_sample_weight_elm():
    """Smoke test - AdaBoostClassifier should work with ELMClassifer."""
    X = Xdigits_binary[:50]
    y = ydigits_binary[:50]

    elm = ELMClassifier(n_hidden=20)
    clf = AdaBoostClassifier(n_estimators=3, base_estimator=elm)
    clf.fit(X, y)
    assert_greater(clf.score(X, y), 0.9)


def test_sparse_matrices():
    """Test that sparse and dense input matrices yield equal output."""
    X = Xdigits_binary[:50]
    y = ydigits_binary[:50]
    X = csr_matrix(X)
    n_hidden = 15
    batch_size = 10

    # Standard ELM
    elm = ELMClassifier(random_state=1, n_hidden=n_hidden)
    # Batch based
    elm_batch_based = ELMClassifier(random_state=1, n_hidden=n_hidden,
                                    batch_size=10)
    # ELM for partial fitting
    elm_parital = ELMClassifier(random_state=1, n_hidden=n_hidden)
    # Train classifiers
    elm.fit(X, y)
    elm_batch_based.fit(X, y)

    for batch_slice in gen_batches(X.shape[0], batch_size):
        elm_parital.partial_fit(X[batch_slice], y[batch_slice])

    # Get decision scores
    y_pred = elm.decision_function(X)
    y_pred_batch_based = elm_batch_based.decision_function(X)
    y_pred_partial = elm_parital.decision_function(X)

    # The prediction values should be the same
    assert_almost_equal(y_pred, y_pred_batch_based)
    assert_almost_equal(y_pred_batch_based, y_pred_partial)


def test_verbose():
    """Test whether verbose works as intended."""
    X = Xboston
    y = yboston

    elm_fit = ELMRegressor(verbose=True)
    elm_batch_fit = ELMRegressor(verbose=True, batch_size=50)
    for elm in [elm_fit, elm_batch_fit]:
        old_stdout = sys.stdout
        sys.stdout = output = StringIO()

        elm.fit(X, y)
        sys.stdout = old_stdout

        assert_not_equal(output.getvalue(), '')


def test_warmstart():
    """Tests that warm_start reuses past solution only if set to true."""
    X = Xboston
    y = yboston

    elm = ELMRegressor(weight_scale=100, warm_start=True)

    elm.fit(X, y)
    first_call_init_weights = elm.coef_hidden_

    elm.weight_scale = 0.001

    elm.fit(X, y)
    second_call_init_weights = elm.coef_hidden_

    # The input-to-hidden weights should stay the same between the two calls.
    assert_array_equal(first_call_init_weights, second_call_init_weights)


def test_weighted_elm():
    """Check class weight impact on AUC.

    Re-weighting classes is expected to improve model performance
    for imbalanced classification problems.
    """
    rng = np.random.RandomState(random_state)
    n_samples_1 = 500
    n_samples_2 = 10
    X = np.r_[1.5 * rng.randn(n_samples_1, 20),
              1.2 * rng.randn(n_samples_2, 20) + [2] * 20]
    y = [0] * (n_samples_1) + [1] * (n_samples_2)

    X_train, X_test, y_train, y_test = cross_validation.train_test_split(
        X, y, test_size=0.8, random_state=random_state)

    n_hidden = 20
    for activation in ACTIVATION_TYPES:
        elm_weightless = ELMClassifier(n_hidden=n_hidden,
                                       class_weight=None,
                                       random_state=random_state)
        elm_weightless.fit(X_train, y_train)

        elm_weight_auto = ELMClassifier(n_hidden=n_hidden,
                                        class_weight='auto',
                                        random_state=random_state)
        elm_weight_auto.fit(X_train, y_train)

        y_pred_weightless = elm_weightless.predict_proba(X_test)[:, 1]
        score_weightless = roc_auc_score(y_test, y_pred_weightless)

        y_pred_weighted = elm_weight_auto.predict_proba(X_test)[:, 1]
        score_weighted = roc_auc_score(y_test, y_pred_weighted)

        assert_greater(score_weighted, score_weightless)
