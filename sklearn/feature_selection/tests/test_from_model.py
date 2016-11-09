import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import skip_if_32bit

from sklearn import datasets
from sklearn.linear_model import LogisticRegression, SGDClassifier, Lasso
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import PassiveAggressiveClassifier

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
    transformer = SelectFromModel(estimator=Lasso(alpha=0.1))
    transformer.fit(X, y)
    X_new = transformer.transform(X)
    mask = np.abs(transformer.estimator_.coef_) > 1e-5
    assert_array_equal(X_new, X[:, mask])


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


def test_calling_fit_reinitializes():
    est = LinearSVC(random_state=0)
    transformer = SelectFromModel(estimator=est)
    transformer.fit(data, y)
    transformer.set_params(estimator__C=100)
    transformer.fit(data, y)
    assert_equal(transformer.estimator_.C, 100)


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
