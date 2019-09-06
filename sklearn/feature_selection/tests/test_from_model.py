import pytest
import numpy as np

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import skip_if_32bit

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


# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
def test_invalid_input():
    clf = SGDClassifier(alpha=0.1, max_iter=10, shuffle=True,
                        random_state=None, tol=None)
    for threshold in ["gobbledigook", ".5 * gobbledigook"]:
        model = SelectFromModel(clf, threshold=threshold)
        model.fit(data, y)
        with pytest.raises(ValueError):
            model.transform(data)


def test_input_estimator_unchanged():
    # Test that SelectFromModel fits on a clone of the estimator.
    est = RandomForestClassifier()
    transformer = SelectFromModel(estimator=est)
    transformer.fit(data, y)
    assert transformer.estimator is est


@pytest.mark.parametrize(
    "max_features, err_type, err_msg",
    [(-1, ValueError, "'max_features' should be 0 and"),
     (data.shape[1] + 1, ValueError, "'max_features' should be 0 and"),
     ('gobbledigook', TypeError, "should be an integer"),
     ('all', TypeError, "should be an integer")]
)
def test_max_features_error(max_features, err_type, err_msg):
    clf = RandomForestClassifier(n_estimators=50, random_state=0)

    transformer = SelectFromModel(estimator=clf,
                                  max_features=max_features,
                                  threshold=-np.inf)
    with pytest.raises(err_type, match=err_msg):
        transformer.fit(data, y)


@pytest.mark.parametrize("max_features", [0, 2, data.shape[1]])
def test_max_features_dim(max_features):
    clf = RandomForestClassifier(n_estimators=50, random_state=0)
    transformer = SelectFromModel(estimator=clf,
                                  max_features=max_features,
                                  threshold=-np.inf)
    X_trans = transformer.fit_transform(data, y)
    assert X_trans.shape[1] == max_features


class FixedImportanceEstimator(BaseEstimator):
    def __init__(self, importances):
        self.importances = importances

    def fit(self, X, y=None):
        self.feature_importances_ = np.array(self.importances)


def test_max_features():
    # Test max_features parameter using various values
    X, y = datasets.make_classification(
        n_samples=1000, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0)
    max_features = X.shape[1]
    est = RandomForestClassifier(n_estimators=50, random_state=0)

    transformer1 = SelectFromModel(estimator=est,
                                   threshold=-np.inf)
    transformer2 = SelectFromModel(estimator=est,
                                   max_features=max_features,
                                   threshold=-np.inf)
    X_new1 = transformer1.fit_transform(X, y)
    X_new2 = transformer2.fit_transform(X, y)
    assert_allclose(X_new1, X_new2)

    # Test max_features against actual model.
    transformer1 = SelectFromModel(estimator=Lasso(alpha=0.025,
                                                   random_state=42))
    X_new1 = transformer1.fit_transform(X, y)
    scores1 = np.abs(transformer1.estimator_.coef_)
    candidate_indices1 = np.argsort(-scores1, kind='mergesort')

    for n_features in range(1, X_new1.shape[1] + 1):
        transformer2 = SelectFromModel(estimator=Lasso(alpha=0.025,
                                       random_state=42),
                                       max_features=n_features,
                                       threshold=-np.inf)
        X_new2 = transformer2.fit_transform(X, y)
        scores2 = np.abs(transformer2.estimator_.coef_)
        candidate_indices2 = np.argsort(-scores2, kind='mergesort')
        assert_allclose(X[:, candidate_indices1[:n_features]],
                        X[:, candidate_indices2[:n_features]])
    assert_allclose(transformer1.estimator_.coef_,
                    transformer2.estimator_.coef_)


def test_max_features_tiebreak():
    # Test if max_features can break tie among feature importance
    X, y = datasets.make_classification(
        n_samples=1000, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0)
    max_features = X.shape[1]

    feature_importances = np.array([4, 4, 4, 4, 3, 3, 3, 2, 2, 1])
    for n_features in range(1, max_features + 1):
        transformer = SelectFromModel(
            FixedImportanceEstimator(feature_importances),
            max_features=n_features,
            threshold=-np.inf)
        X_new = transformer.fit_transform(X, y)
        selected_feature_indices = np.where(transformer._get_support_mask())[0]
        assert_array_equal(selected_feature_indices, np.arange(n_features))
        assert X_new.shape[1] == n_features


def test_threshold_and_max_features():
    X, y = datasets.make_classification(
        n_samples=1000, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0)
    est = RandomForestClassifier(n_estimators=50, random_state=0)

    transformer1 = SelectFromModel(estimator=est, max_features=3,
                                   threshold=-np.inf)
    X_new1 = transformer1.fit_transform(X, y)

    transformer2 = SelectFromModel(estimator=est, threshold=0.04)
    X_new2 = transformer2.fit_transform(X, y)

    transformer3 = SelectFromModel(estimator=est, max_features=3,
                                   threshold=0.04)
    X_new3 = transformer3.fit_transform(X, y)
    assert X_new3.shape[1] == min(X_new1.shape[1], X_new2.shape[1])
    selected_indices = transformer3.transform(
        np.arange(X.shape[1])[np.newaxis, :])
    assert_allclose(X_new3, X[:, selected_indices[0]])


@skip_if_32bit
def test_feature_importances():
    X, y = datasets.make_classification(
        n_samples=1000, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0)

    est = RandomForestClassifier(n_estimators=50, random_state=0)
    for threshold, func in zip(["mean", "median"], [np.mean, np.median]):
        transformer = SelectFromModel(estimator=est, threshold=threshold)
        transformer.fit(X, y)
        assert hasattr(transformer.estimator_, 'feature_importances_')

        X_new = transformer.transform(X)
        assert X_new.shape[1] < X.shape[1]
        importances = transformer.estimator_.feature_importances_

        feature_mask = np.abs(importances) > func(importances)
        assert_array_almost_equal(X_new, X[:, feature_mask])


def test_sample_weight():
    # Ensure sample weights are passed to underlying estimator
    X, y = datasets.make_classification(
        n_samples=100, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0)

    # Check with sample weights
    sample_weight = np.ones(y.shape)
    sample_weight[y == 1] *= 100

    est = LogisticRegression(random_state=0, fit_intercept=False)
    transformer = SelectFromModel(estimator=est)
    transformer.fit(X, y, sample_weight=None)
    mask = transformer._get_support_mask()
    transformer.fit(X, y, sample_weight=sample_weight)
    weighted_mask = transformer._get_support_mask()
    assert not np.all(weighted_mask == mask)
    transformer.fit(X, y, sample_weight=3 * sample_weight)
    reweighted_mask = transformer._get_support_mask()
    assert np.all(weighted_mask == reweighted_mask)


def test_coef_default_threshold():
    X, y = datasets.make_classification(
        n_samples=100, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0)

    # For the Lasso and related models, the threshold defaults to 1e-5
    transformer = SelectFromModel(estimator=Lasso(alpha=0.1,
                                  random_state=42))
    transformer.fit(X, y)
    X_new = transformer.transform(X)
    mask = np.abs(transformer.estimator_.coef_) > 1e-5
    assert_array_almost_equal(X_new, X[:, mask])


@skip_if_32bit
def test_2d_coef():
    X, y = datasets.make_classification(
        n_samples=1000, n_features=10, n_informative=3, n_redundant=0,
        n_repeated=0, shuffle=False, random_state=0, n_classes=4)

    est = LogisticRegression()
    for threshold, func in zip(["mean", "median"], [np.mean, np.median]):
        for order in [1, 2, np.inf]:
            # Fit SelectFromModel a multi-class problem
            transformer = SelectFromModel(estimator=LogisticRegression(),
                                          threshold=threshold,
                                          norm_order=order)
            transformer.fit(X, y)
            assert hasattr(transformer.estimator_, 'coef_')
            X_new = transformer.transform(X)
            assert X_new.shape[1] < X.shape[1]

            # Manually check that the norm is correctly performed
            est.fit(X, y)
            importances = np.linalg.norm(est.coef_, axis=0, ord=order)
            feature_mask = importances > func(importances)
            assert_array_almost_equal(X_new, X[:, feature_mask])


# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
def test_partial_fit():
    est = PassiveAggressiveClassifier(random_state=0, shuffle=False,
                                      max_iter=5, tol=None)
    transformer = SelectFromModel(estimator=est)
    transformer.partial_fit(data, y,
                            classes=np.unique(y))
    old_model = transformer.estimator_
    transformer.partial_fit(data, y,
                            classes=np.unique(y))
    new_model = transformer.estimator_
    assert old_model is new_model

    X_transform = transformer.transform(data)
    transformer.fit(np.vstack((data, data)), np.concatenate((y, y)))
    assert_array_almost_equal(X_transform, transformer.transform(data))

    # check that if est doesn't have partial_fit, neither does SelectFromModel
    transformer = SelectFromModel(estimator=RandomForestClassifier())
    assert not hasattr(transformer, "partial_fit")


def test_calling_fit_reinitializes():
    est = LinearSVC(random_state=0)
    transformer = SelectFromModel(estimator=est)
    transformer.fit(data, y)
    transformer.set_params(estimator__C=100)
    transformer.fit(data, y)
    assert transformer.estimator_.C == 100


# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
def test_prefit():
    # Test all possible combinations of the prefit parameter.

    # Passing a prefit parameter with the selected model
    # and fitting a unfit model with prefit=False should give same results.
    clf = SGDClassifier(alpha=0.1, max_iter=10, shuffle=True,
                        random_state=0, tol=None)
    model = SelectFromModel(clf)
    model.fit(data, y)
    X_transform = model.transform(data)
    clf.fit(data, y)
    model = SelectFromModel(clf, prefit=True)
    assert_array_almost_equal(model.transform(data), X_transform)

    # Check that the model is rewritten if prefit=False and a fitted model is
    # passed
    model = SelectFromModel(clf, prefit=False)
    model.fit(data, y)
    assert_array_almost_equal(model.transform(data), X_transform)

    # Check that prefit=True and calling fit raises a ValueError
    model = SelectFromModel(clf, prefit=True)
    with pytest.raises(ValueError):
        model.fit(data, y)


def test_threshold_string():
    est = RandomForestClassifier(n_estimators=50, random_state=0)
    model = SelectFromModel(est, threshold="0.5*mean")
    model.fit(data, y)
    X_transform = model.transform(data)

    # Calculate the threshold from the estimator directly.
    est.fit(data, y)
    threshold = 0.5 * np.mean(est.feature_importances_)
    mask = est.feature_importances_ > threshold
    assert_array_almost_equal(X_transform, data[:, mask])


# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
def test_threshold_without_refitting():
    # Test that the threshold can be set without refitting the model.
    clf = SGDClassifier(alpha=0.1, max_iter=10, shuffle=True,
                        random_state=0, tol=None)
    model = SelectFromModel(clf, threshold="0.1 * mean")
    model.fit(data, y)
    X_transform = model.transform(data)

    # Set a higher threshold to filter out more features.
    model.threshold = "1.0 * mean"
    assert X_transform.shape[1] > model.transform(data).shape[1]
