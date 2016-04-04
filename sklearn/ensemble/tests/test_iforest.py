"""
Testing for Isolation Forest algorithm (sklearn.ensemble.iforest).
"""

# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import ignore_warnings

from sklearn.grid_search import ParameterGrid
from sklearn.ensemble import IsolationForest
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_boston, load_iris
from sklearn.utils import check_random_state
from sklearn.metrics import roc_auc_score

from scipy.sparse import csc_matrix, csr_matrix

rng = check_random_state(0)

# load the iris dataset
# and randomly permute it
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_iforest():
    """Check Isolation Forest for various parameter settings."""
    X_train = np.array([[0, 1], [1, 2]])
    X_test = np.array([[2, 1], [1, 1]])

    grid = ParameterGrid({"n_estimators": [3],
                          "max_samples": [0.5, 1.0, 3],
                          "bootstrap": [True, False]})

    with ignore_warnings():
        for params in grid:
            IsolationForest(random_state=rng,
                            **params).fit(X_train).predict(X_test)


def test_iforest_sparse():
    """Check IForest for various parameter settings on sparse input."""
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data[:50],
                                                        boston.target[:50],
                                                        random_state=rng)
    grid = ParameterGrid({"max_samples": [0.5, 1.0],
                          "bootstrap": [True, False]})

    for sparse_format in [csc_matrix, csr_matrix]:
        X_train_sparse = sparse_format(X_train)
        X_test_sparse = sparse_format(X_test)

        for params in grid:
            # Trained on sparse format
            sparse_classifier = IsolationForest(
                n_estimators=10, random_state=1, **params).fit(X_train_sparse)
            sparse_results = sparse_classifier.predict(X_test_sparse)

            # Trained on dense format
            dense_classifier = IsolationForest(
                n_estimators=10, random_state=1, **params).fit(X_train)
            dense_results = dense_classifier.predict(X_test)

            assert_array_equal(sparse_results, dense_results)
            assert_array_equal(sparse_results, dense_results)


def test_iforest_error():
    """Test that it gives proper exception on deficient input."""
    X = iris.data

    # Test max_samples
    assert_raises(ValueError,
                  IsolationForest(max_samples=-1).fit, X)
    assert_raises(ValueError,
                  IsolationForest(max_samples=0.0).fit, X)
    assert_raises(ValueError,
                  IsolationForest(max_samples=2.0).fit, X)
    # The dataset has less than 256 samples, explicitly setting max_samples > n_samples
    # should result in a warning. If not set explicitly there should be no warning
    assert_warns_message(UserWarning,
                         "max_samples will be set to n_samples for estimation",
                         IsolationForest(max_samples=1000).fit, X)
    assert_no_warnings(IsolationForest(max_samples='auto').fit, X)
    assert_raises(ValueError,
                  IsolationForest(max_samples='foobar').fit, X)


def test_recalculate_max_depth():
    """Check that max_depth is recalculated when max_samples is reset to n_samples"""
    X = iris.data
    clf = IsolationForest().fit(X)
    for est in clf.estimators_:
        assert_equal(est.max_depth, int(np.ceil(np.log2(X.shape[0]))))


def test_max_samples_attribute():
    X = iris.data
    clf = IsolationForest().fit(X)
    assert_equal(clf.max_samples_, X.shape[0])

    clf = IsolationForest(max_samples=500)
    assert_warns_message(UserWarning,
                         "max_samples will be set to n_samples for estimation",
                         clf.fit, X)
    assert_equal(clf.max_samples_, X.shape[0])

    clf = IsolationForest(max_samples=0.4).fit(X)
    assert_equal(clf.max_samples_, 0.4*X.shape[0])


def test_iforest_parallel_regression():
    """Check parallel regression."""
    rng = check_random_state(0)

    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    ensemble = IsolationForest(n_jobs=3,
                               random_state=0).fit(X_train)

    ensemble.set_params(n_jobs=1)
    y1 = ensemble.predict(X_test)
    ensemble.set_params(n_jobs=2)
    y2 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y2)

    ensemble = IsolationForest(n_jobs=1,
                               random_state=0).fit(X_train)

    y3 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y3)


def test_iforest_performance():
    """Test Isolation Forest performs well"""

    # Generate train/test data
    rng = check_random_state(2)
    X = 0.3 * rng.randn(120, 2)
    X_train = np.r_[X + 2, X - 2]
    X_train = X[:100]

    # Generate some abnormal novel observations
    X_outliers = rng.uniform(low=-4, high=4, size=(20, 2))
    X_test = np.r_[X[100:], X_outliers]
    y_test = np.array([0] * 20 + [1] * 20)

    # fit the model
    clf = IsolationForest(max_samples=100, random_state=rng).fit(X_train)

    # predict scores (the lower, the more normal)
    y_pred = clf.predict(X_test)

    # check that there is at most 6 errors (false positive or false negative)
    assert_greater(roc_auc_score(y_test, y_pred), 0.98)


def test_iforest_works():
    # toy sample (the last two samples are outliers)
    X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1], [6, 3], [-4, 7]]

    # Test LOF
    clf = IsolationForest(random_state=rng)
    clf.fit(X)
    pred = clf.predict(X)

    # assert detect outliers:
    assert_greater(np.min(pred[-2:]), np.max(pred[:-2]))
