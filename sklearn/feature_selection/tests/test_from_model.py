import numpy as np
import scipy.sparse as sp

from nose.tools import assert_raises, assert_true

from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import skip_if_32bit
from sklearn.utils.testing import assert_equal

from sklearn import datasets
from sklearn.linear_model import LogisticRegression, SGDClassifier, Lasso
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.base import BaseEstimator

iris = datasets.load_iris()
data, y = iris.data, iris.target
rng = np.random.RandomState(0)


def test_transform_linear_model():
    for clf in (LogisticRegression(C=0.1),
                LinearSVC(C=0.01, dual=False),
                SGDClassifier(alpha=0.001, n_iter=50, shuffle=True,
                              random_state=0)):
        for thresh in (None, ".09*mean", "1e-5 * median"):
            for func in (np.array, sp.csr_matrix):
                X = func(data)
                clf.set_params(penalty="l1")
                clf.fit(X, y)
                X_new = assert_warns(
                    DeprecationWarning, clf.transform, X, thresh)
                if isinstance(clf, SGDClassifier):
                    assert_true(X_new.shape[1] <= X.shape[1])
                else:
                    assert_less(X_new.shape[1], X.shape[1])
                clf.set_params(penalty="l2")
                clf.fit(X_new, y)
                pred = clf.predict(X_new)
                assert_greater(np.mean(pred == y), 0.7)


def test_invalid_input():
    clf = SGDClassifier(alpha=0.1, n_iter=10, shuffle=True, random_state=None)
    for threshold in ["gobbledigook", ".5 * gobbledigook"]:
        model = SelectFromModel(clf, threshold=threshold)
        model.fit(data, y)
        assert_raises(ValueError, model.transform, data)


def test_input_estimator_unchanged():
    """
    Test that SelectFromModel fits on a clone of the estimator.
    """
    est = RandomForestClassifier()
    transformer = SelectFromModel(estimator=est)
    transformer.fit(data, y)
    assert_true(transformer.estimator is est)


def check_invalid_max_features(est, X, y):
    max_features = X.shape[1]
    for invalid_max_n_feature in [-1, max_features + 1, 'gobbledigook']:
        transformer = SelectFromModel(estimator=est,
                                      max_features=invalid_max_n_feature)
        assert_raises(ValueError, transformer.fit, X, y)


def check_valid_max_features(est, X, y):
    max_features = X.shape[1]
    for valid_max_n_feature in [0, max_features, 'all', 5]:
        transformer = SelectFromModel(estimator=est,
                                      max_features=valid_max_n_feature)
        X_new = transformer.fit_transform(X, y)
        if valid_max_n_feature == 'all':
            valid_max_n_feature = max_features
        assert_equal(X_new.shape[1], valid_max_n_feature)


class FixedImportanceEstimator(BaseEstimator):
    def __init__(self, importances):
        self.importances = importances

    def fit(self, X, y=None):
        self.feature_importances_ = np.array(self.importances)


def check_max_features(est, X, y):
    X = X.copy()
    max_features = X.shape[1]

    check_valid_max_features(est, X, y)
    check_invalid_max_features(est, X, y)

    transformer1 = SelectFromModel(estimator=est, max_features='all')
    transformer2 = SelectFromModel(estimator=est,
                                   max_features=max_features)
    X_new1 = transformer1.fit_transform(X, y)
    X_new2 = transformer2.fit_transform(X, y)
    assert_array_equal(X_new1, X_new2)

    # Test max_features against actual model.

    transformer1 = SelectFromModel(estimator=Lasso(alpha=0.025))
    X_new1 = transformer1.fit_transform(X, y)
    for n_features in range(1, X_new1.shape[1] + 1):
        transformer2 = SelectFromModel(estimator=Lasso(alpha=0.025),
                                       max_features=n_features)
        X_new2 = transformer2.fit_transform(X, y)
        assert_array_equal(X_new1[:, :n_features], X_new2)
        assert_array_equal(transformer1.estimator_.coef_,
                           transformer2.estimator_.coef_)

    # Test if max_features can break tie among feature importance

    feature_importances = np.array([4, 4, 4, 4, 3, 3, 3, 2, 2, 1])
    for n_features in range(1, max_features + 1):
        transformer = SelectFromModel(
            FixedImportanceEstimator(feature_importances),
            max_features=n_features)
        X_new = transformer.fit_transform(X, y)
        selected_feature_indices = np.where(transformer._get_support_mask())[0]
        assert_array_equal(selected_feature_indices, np.arange(n_features))
        assert_equal(X_new.shape[1], n_features)


def check_threshold_and_max_features(est, X, y):
    transformer1 = SelectFromModel(estimator=est, max_features=3)
    X_new1 = transformer1.fit_transform(X, y)

    transformer2 = SelectFromModel(estimator=est, threshold=0.04)
    X_new2 = transformer2.fit_transform(X, y)

    transformer3 = SelectFromModel(estimator=est, max_features=3,
                                   threshold=0.04)
    X_new3 = transformer3.fit_transform(X, y)
    assert_equal(X_new3.shape[1], min(X_new1.shape[1], X_new2.shape[1]))
    selected_indices = \
        transformer3.transform(np.arange(X.shape[1]))[np.newaxis, :]
    assert_array_equal(X_new3, X[:, selected_indices[0][0]])

    """
    If threshold and max_features are not provided, all features are
    returned, use threshold=None if it is not required.
    """
    transformer = SelectFromModel(estimator=Lasso(alpha=0.1))
    X_new = transformer.fit_transform(X, y)
    assert_array_equal(X, X_new)

    transformer = SelectFromModel(estimator=Lasso(alpha=0.1), max_features=3)
    X_new = transformer.fit_transform(X, y)
    assert_equal(X_new.shape[1], 3)

    # Threshold will be applied if it is not None
    transformer = SelectFromModel(estimator=Lasso(alpha=0.1), threshold=1e-5)
    X_new = transformer.fit_transform(X, y)
    mask = np.abs(transformer.estimator_.coef_) > 1e-5
    assert_array_equal(X_new, X[:, mask])

    transformer = SelectFromModel(estimator=Lasso(alpha=0.1), threshold=1e-5,
                                  max_features=4)
    X_new = transformer.fit_transform(X, y)
    mask = np.abs(transformer.estimator_.coef_) > 1e-5
    assert_array_equal(X_new, X[:, mask])


@skip_if_32bit
def test_feature_importances():
    X, y = datasets.make_classification(
        n_samples=1000, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0)

    est = RandomForestClassifier(n_estimators=50, random_state=0)
    for threshold, func in zip(["mean", "median"], [np.mean, np.median]):
        transformer = SelectFromModel(estimator=est, threshold=threshold)
        transformer.fit(X, y)
        assert_true(hasattr(transformer.estimator_, 'feature_importances_'))

        X_new = transformer.transform(X)
        assert_less(X_new.shape[1], X.shape[1])
        importances = transformer.estimator_.feature_importances_

        feature_mask = np.abs(importances) > func(importances)
        assert_array_almost_equal(X_new, X[:, feature_mask])

    # Check with sample weights
    sample_weight = np.ones(y.shape)
    sample_weight[y == 1] *= 100

    est = RandomForestClassifier(n_estimators=50, random_state=0)
    transformer = SelectFromModel(estimator=est)
    transformer.fit(X, y, sample_weight=sample_weight)
    importances = transformer.estimator_.feature_importances_
    transformer.fit(X, y, sample_weight=3 * sample_weight)
    importances_bis = transformer.estimator_.feature_importances_
    assert_almost_equal(importances, importances_bis)

    # For the Lasso and related models, the threshold defaults to 1e-5
    transformer = SelectFromModel(estimator=Lasso(alpha=0.1), threshold=1e-5)
    transformer.fit(X, y)
    X_new = transformer.transform(X)
    mask = np.abs(transformer.estimator_.coef_) > 1e-5
    assert_array_equal(X_new, X[:, mask])

    # Test max_features parameter using various values
    check_max_features(est, X, y)
    check_threshold_and_max_features(est, X, y)


def test_partial_fit():
    est = PassiveAggressiveClassifier(random_state=0, shuffle=False)
    transformer = SelectFromModel(estimator=est)
    transformer.partial_fit(data, y,
                            classes=np.unique(y))
    old_model = transformer.estimator_
    transformer.partial_fit(data, y,
                            classes=np.unique(y))
    new_model = transformer.estimator_
    assert_true(old_model is new_model)

    X_transform = transformer.transform(data)
    transformer.fit(np.vstack((data, data)), np.concatenate((y, y)))
    assert_array_equal(X_transform, transformer.transform(data))


def test_warm_start():
    est = PassiveAggressiveClassifier(warm_start=True, random_state=0)
    transformer = SelectFromModel(estimator=est)
    transformer.fit(data, y)
    old_model = transformer.estimator_
    transformer.fit(data, y)
    new_model = transformer.estimator_
    assert_true(old_model is new_model)


def test_prefit():
    """
    Test all possible combinations of the prefit parameter.
    """
    # Passing a prefit parameter with the selected model
    # and fitting a unfit model with prefit=False should give same results.
    clf = SGDClassifier(alpha=0.1, n_iter=10, shuffle=True, random_state=0)
    model = SelectFromModel(clf)
    model.fit(data, y)
    X_transform = model.transform(data)
    clf.fit(data, y)
    model = SelectFromModel(clf, prefit=True)
    assert_array_equal(model.transform(data), X_transform)

    # Check that the model is rewritten if prefit=False and a fitted model is
    # passed
    model = SelectFromModel(clf, prefit=False)
    model.fit(data, y)
    assert_array_equal(model.transform(data), X_transform)

    # Check that prefit=True and calling fit raises a ValueError
    model = SelectFromModel(clf, prefit=True)
    assert_raises(ValueError, model.fit, data, y)


def test_threshold_string():
    est = RandomForestClassifier(n_estimators=50, random_state=0)
    model = SelectFromModel(est, threshold="0.5*mean")
    model.fit(data, y)
    X_transform = model.transform(data)

    # Calculate the threshold from the estimator directly.
    est.fit(data, y)
    threshold = 0.5 * np.mean(est.feature_importances_)
    mask = est.feature_importances_ > threshold
    assert_array_equal(X_transform, data[:, mask])


def test_threshold_without_refitting():
    """Test that the threshold can be set without refitting the model."""
    clf = SGDClassifier(alpha=0.1, n_iter=10, shuffle=True, random_state=0)
    model = SelectFromModel(clf, threshold=0.1)
    model.fit(data, y)
    X_transform = model.transform(data)

    # Set a higher threshold to filter out more features.
    model.threshold = 1.0
    assert_greater(X_transform.shape[1], model.transform(data).shape[1])
