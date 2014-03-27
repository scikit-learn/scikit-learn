"""
Testing for the SelectFromModel (sklearn.preprocessing.select_from_model).
"""

# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>
# License: BSD 3 clause

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true

from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import AdaBoostClassifier
from sklearn.dummy import DummyClassifier
from sklearn import datasets

iris = datasets.load_iris()
X = iris.data
y = iris.target


def test_validate_estimator():
    # Test `None` as estimator
    assert_raises(ValueError, SelectFromModel, estimator=None)
    assert_raises(ValueError, SelectFromModel)

    est = AdaBoostClassifier()
    tfm = SelectFromModel(estimator=est)
    assert_equal(tfm.estimator_, est)


def test_fit():
    est = DummyClassifier()
    tfm = SelectFromModel(estimator=est)
    assert_raises(ValueError, tfm.fit, X, y)


def test_feature_importances():
    est = AdaBoostClassifier()
    tfm = SelectFromModel(estimator=est)

    tfm.fit(X, y)
    assert_true(hasattr(tfm, 'feature_importances_'))

    X_new = tfm.transform(X, threshold="mean")
    assert_less(X_new.shape[1], X.shape[1])

    feature_mask = tfm.feature_importances_ > tfm.feature_importances_.mean(
    )
    assert_array_almost_equal(X_new, X[:, feature_mask])
