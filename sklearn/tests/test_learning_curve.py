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
    def __init__(self, n_max_train_samples):
        self.n_max_train_samples = n_max_train_samples
        self.n_train_samples = 0
        self.X_subset = None

    def fit(self, X_subset, y_subset):
        self.X_subset = X_subset
        self.n_train_samples = X_subset.shape[0]
        return self

    def predict(self, X):
        raise NotImplementedError

    def score(self, X=None, Y=None):
        # training score becomes worse (2 -> 1), test error better (0 -> 1)
        if self._is_training_data(X):
            return 2. - float(self.n_train_samples) / self.n_max_train_samples
        else:
            return float(self.n_train_samples) / self.n_max_train_samples

    def _is_training_data(self, X):
        return X is self.X_subset

    def get_params(self, deep=False):
        return {"n_max_train_samples" : self.n_max_train_samples}

    def set_params(self, **params):
        self.n_max_train_samples = params["n_max_train_samples"]
        return self


class MockIncrementalImprovingClassifier(MockImprovingClassifier):
    """Dummy classifier that provides partial_fit"""
    def __init__(self, n_max_train_samples):
        super(MockIncrementalImprovingClassifier, self).__init__(
                n_max_train_samples)
        self.x = None

    def _is_training_data(self, X):
        return self.x in X

    def partial_fit(self, X, y, **params):
        self.n_train_samples += X.shape[0]
        self.x = X[0]


def test_learning_curve():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(20)
    n_samples_range, train_scores, test_scores = learning_curve(estimator,
                                                                X, y, cv=3)
    assert_array_equal(n_samples_range, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores, np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores, np.linspace(0.1, 1.0, 10))


def test_incremental_learning_not_possible():
    X, y = make_classification(n_samples=2, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    # The mockup does not have partial_fit()
    estimator = MockImprovingClassifier(1)
    assert_raises(ValueError, learning_curve, estimator, X, y,
                  exploit_incremental_learning=True)


def test_incremental_learning():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockIncrementalImprovingClassifier(20)
    n_samples_range, train_scores, test_scores = learning_curve(
            estimator, X, y, cv=3, exploit_incremental_learning=True)
    assert_array_equal(n_samples_range, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores, np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores, np.linspace(0.1, 1.0, 10))


def test_batch_and_incremental_learning_are_equal():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = PassiveAggressiveClassifier(n_iter=1, shuffle=False)
    n_samples_range_inc, train_scores_inc, test_scores_inc = learning_curve(
            estimator, X, y, cv=3, exploit_incremental_learning=True)
    n_samples_range_batch, train_scores_batch, test_scores_batch = \
            learning_curve(estimator, X, y, cv=3,
                           exploit_incremental_learning=False)
    assert_array_equal(n_samples_range_inc, n_samples_range_batch)
    assert_array_almost_equal(train_scores_inc, train_scores_batch)
    assert_array_almost_equal(test_scores_inc, test_scores_batch)


def test_n_sample_range_out_of_bounds():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(20)
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  n_samples_range=[0.0, 1.0])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  n_samples_range=[0.1, 1.1])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  n_samples_range=[0, 20])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  n_samples_range=[1, 21])


def test_remove_duplicate_sample_sizes():
    X, y = make_classification(n_samples=3, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(2)
    n_samples_range, _, _ = assert_warns(RuntimeWarning,
            learning_curve, estimator, X, y, cv=3,
            n_samples_range=np.linspace(0.33, 1.0, 3))
    assert_array_equal(n_samples_range, [1, 2])


def test_learning_curve_with_boolean_indices():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(20)
    cv = KFold(n=30, n_folds=3, indices=False)
    n_samples_range, train_scores, test_scores = learning_curve(estimator,
                                                                X, y, cv=cv)
    assert_array_equal(n_samples_range, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores, np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores, np.linspace(0.1, 1.0, 10))
    
