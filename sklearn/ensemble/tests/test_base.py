"""
Testing for the base module (sklearn.ensemble.base).
"""

# Authors: Gilles Louppe
# License: BSD 3

from numpy.testing import assert_equal
from nose.tools import assert_raises

from sklearn.ensemble import BaseEnsemble
from sklearn.tree import DecisionTreeClassifier


def test_base():
    """Check BaseEnsemble methods."""
    tree = DecisionTreeClassifier()
    ensemble = BaseEnsemble(base_estimator=tree, n_estimators=3)

    ensemble._make_estimator()
    ensemble._make_estimator()
    ensemble._make_estimator()
    ensemble._make_estimator(append=False)

    assert_equal(3, len(ensemble))
    assert_equal(3, len(ensemble.estimators_))

    assert isinstance(ensemble[0], DecisionTreeClassifier)


def test_error():
    """Check that proper errors are triggered."""
    def instantiate(class_name, **params):
        return class_name(**params)

    base_estimator = object()
    assert_raises(TypeError, instantiate, class_name=BaseEnsemble,
                  base_estimator=base_estimator, n_estimators=1)

    base_estimator = DecisionTreeClassifier()
    assert_raises(ValueError, instantiate, class_name=BaseEnsemble,
                  base_estimator=base_estimator, n_estimators=-1)
