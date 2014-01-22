# Author: Alexander Fabisch <afabisch@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import sys
from sklearn.externals.six.moves import cStringIO as StringIO
import numpy as np
from sklearn.learning_curve import learning_curve
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.datasets import make_classification
from sklearn.cross_validation import KFold
from sklearn.linear_model import PassiveAggressiveClassifier


class MockImprovingClassifier(object):
    """Dummy classifier to test the learning curve"""
    def __init__(self, n_max_train_sizes):
        self.n_max_train_sizes = n_max_train_sizes
        self.train_sizes = 0
        self.X_subset = None

    def fit(self, X_subset, y_subset):
        self.X_subset = X_subset
        self.train_sizes = X_subset.shape[0]
        return self

    def predict(self, X):
        raise NotImplementedError

    def score(self, X=None, Y=None):
        # training score becomes worse (2 -> 1), test error better (0 -> 1)
        if self._is_training_data(X):
            return 2. - float(self.train_sizes) / self.n_max_train_sizes
        else:
            return float(self.train_sizes) / self.n_max_train_sizes

    def _is_training_data(self, X):
        return X is self.X_subset

    def get_params(self, deep=False):
        return {"n_max_train_sizes": self.n_max_train_sizes}

    def set_params(self, **params):
        self.n_max_train_sizes = params["n_max_train_sizes"]
        return self


class MockIncrementalImprovingClassifier(MockImprovingClassifier):
    """Dummy classifier that provides partial_fit"""
    def __init__(self, n_max_train_sizes):
        super(MockIncrementalImprovingClassifier, self).__init__(
            n_max_train_sizes)
        self.x = None

    def _is_training_data(self, X):
        return self.x in X

    def partial_fit(self, X, y, **params):
        self.train_sizes += X.shape[0]
        self.x = X[0]


def test_learning_curve():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(20)
    train_sizes, train_scores, test_scores = learning_curve(estimator, X, y,
                                                              cv=3)
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores, np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores, np.linspace(0.1, 1.0, 10))


def test_learning_curve_verbose():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(20)

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        train_sizes, train_scores, test_scores = \
            learning_curve(estimator, X, y, cv=3, verbose=1)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    assert("[learning_curve]" in out)


def test_learning_curve_incremental_learning_not_possible():
    X, y = make_classification(n_samples=2, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    # The mockup does not have partial_fit()
    estimator = MockImprovingClassifier(1)
    assert_raises(ValueError, learning_curve, estimator, X, y,
                  exploit_incremental_learning=True)


def test_learning_curve_incremental_learning():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockIncrementalImprovingClassifier(20)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=3, exploit_incremental_learning=True)
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores, np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores, np.linspace(0.1, 1.0, 10))


def test_learning_curve_batch_and_incremental_learning_are_equal():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    train_sizes = np.linspace(0.2, 1.0, 5)
    estimator = PassiveAggressiveClassifier(n_iter=1, shuffle=False)

    train_sizes_inc, train_scores_inc, test_scores_inc = \
        learning_curve(
            estimator, X, y, train_sizes=train_sizes,
            cv=3, exploit_incremental_learning=True)
    train_sizes_batch, train_scores_batch, test_scores_batch = \
        learning_curve(
            estimator, X, y, cv=3, train_sizes=train_sizes,
            exploit_incremental_learning=False)

    assert_array_equal(train_sizes_inc, train_sizes_batch)
    assert_array_almost_equal(train_scores_inc, train_scores_batch)
    assert_array_almost_equal(test_scores_inc, test_scores_batch)


def test_learning_curve_n_sample_range_out_of_bounds():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(20)
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0.0, 1.0])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0.1, 1.1])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0, 20])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[1, 21])


def test_learning_curve_remove_duplicate_sample_sizes():
    X, y = make_classification(n_samples=3, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(2)
    train_sizes, _, _ = assert_warns(
        RuntimeWarning, learning_curve, estimator, X, y, cv=3,
        train_sizes=np.linspace(0.33, 1.0, 3))
    assert_array_equal(train_sizes, [1, 2])


def test_learning_curve_with_boolean_indices():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(20)
    cv = KFold(n=30, n_folds=3, indices=False)
    train_sizes, train_scores, test_scores = learning_curve(estimator, X, y,
                                                            cv=cv)
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores, np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores, np.linspace(0.1, 1.0, 10))
