import numpy as np
from sklearn.learning_curve import learning_curve
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.datasets import make_classification
from sklearn.svm import SVC

class MockImprovingClassifier(object):
    """Dummy classifier to test the learning curve"""
    def __init__(self, n_max_train_samples):
        self.n_max_train_samples = n_max_train_samples
        self.n_train_samples = 0

    def fit(self, X_subset, y_subset):
        self.X_subset = X_subset
        self.y_subset = y_subset
        self.n_train_samples = X_subset.shape[0]
        return self

    def predict(self, X):
        raise NotImplementedError

    def score(self, X=None, Y=None):
        # training score becomes worse (2 -> 1), test error better (0 -> 1)
        if X is self.X_subset:
            return 2. - float(self.n_train_samples) / self.n_max_train_samples
        else:
            return float(self.n_train_samples) / self.n_max_train_samples

    def get_params(self, deep=False):
        return {"n_max_train_samples" : self.n_max_train_samples}

    def set_params(self, **params):
        self.n_max_train_samples = params["n_max_train_samples"]
        return self


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

def test_remove_multiple_sample_sizes():
    X, y = make_classification(n_samples=3, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingClassifier(2)
    n_samples_range, _, _ = assert_warns(RuntimeWarning,
            learning_curve, estimator, X, y, cv=3,
            n_samples_range=np.linspace(0.33, 1.0, 3))
    assert_array_equal(n_samples_range, [1, 2])
