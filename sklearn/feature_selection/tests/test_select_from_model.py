"""
Testing for the SelectFromModel (sklearn.feature_selection.SelectFromModel).
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
from sklearn import datasets

iris = datasets.load_iris()
X = iris.data
y = iris.target


def test_validate_estimator():
    #Test `None` as estimator
    transformer = SelectFromModel(estimator=None)
    assert_raises(ValueError, transformer.fit, X, y)

    est = AdaBoostClassifier()
    transformer = SelectFromModel(estimator=est)
    transformer.fit(X, y)
    assert_equal(transformer.estimator_, est)


def test_feature_importances():
    est = AdaBoostClassifier()
    transformer = SelectFromModel(estimator=est)

    transformer.fit(X, y)
    assert_true(hasattr(transformer.estimator, 'feature_importances_'))

    X_new = transformer.transform(X, threshold="mean")
    assert_less(X_new.shape[1], X.shape[1])

    feature_mask = (transformer.estimator.feature_importances_ >
                    transformer.estimator.feature_importances_.mean())
    assert_array_almost_equal(X_new, X[:, feature_mask])
