"""
Testing for the base module (sklearn.ensemble.base).
"""

import numpy as np
from numpy.testing import assert_equal
from nose.tools import assert_raises

from sklearn.ensemble import BaseEnsemble
from sklearn.tree import DecisionTreeClassifier


def test_base():
    """Check BaseEnsemble methods."""
    tree = DecisionTreeClassifier()
    ensemble = BaseEnsemble(base_estimator=tree, n_estimators=1)
    ensemble.estimators.append(tree)

    assert_equal(1, len(ensemble))
    assert_equal(tree, ensemble[0])


def test_error():
    """Check that proper errors are triggered."""
    def instantiate(class_name, **params):
        return class_name(**params)

    base_estimator = object()
    assert_raises(TypeError, instantiate, class_name=BaseEnsemble, base_estimator=base_estimator, n_estimators=1)

    base_estimator = DecisionTreeClassifier()
    assert_raises(ValueError, instantiate, class_name=BaseEnsemble, base_estimator=base_estimator, n_estimators=-1)

