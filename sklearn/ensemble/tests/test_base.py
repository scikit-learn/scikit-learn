"""
Testing for the base module (sklearn.ensemble.base).
"""

# Authors: Gilles Louppe
# License: BSD 3 clause

import numpy as np
from numpy.testing import assert_equal
from nose.tools import assert_true

from sklearn.utils.testing import assert_raise_message
from sklearn.datasets import load_iris
from sklearn.ensemble import BaggingClassifier
from sklearn.linear_model import Perceptron


def test_base():
    # Check BaseEnsemble methods.
    ensemble = BaggingClassifier(base_estimator=Perceptron(), n_estimators=3)

    iris = load_iris()
    ensemble.fit(iris.data, iris.target)
    ensemble.estimators_ = []  # empty the list and create estimators manually

    ensemble._make_estimator()
    ensemble._make_estimator()
    ensemble._make_estimator()
    ensemble._make_estimator(append=False)

    assert_equal(3, len(ensemble))
    assert_equal(3, len(ensemble.estimators_))

    assert_true(isinstance(ensemble[0], Perceptron))

    np_int_ensemble = BaggingClassifier(base_estimator=Perceptron(),
                                        n_estimators=np.int32(3))
    np_int_ensemble.fit(iris.data, iris.target)


def test_base_zero_n_estimators():
    # Check that instantiating a BaseEnsemble with n_estimators<=0 raises
    # a ValueError.
    ensemble = BaggingClassifier(base_estimator=Perceptron(),
                                 n_estimators=0)
    iris = load_iris()
    assert_raise_message(ValueError,
                         "n_estimators must be greater than zero, got 0.",
                         ensemble.fit, iris.data, iris.target)


def test_base_not_int_n_estimators():
    # Check that instantiating a BaseEnsemble with a string as n_estimators
    # raises a ValueError demanding n_estimators to be supplied as an integer.
    string_ensemble = BaggingClassifier(base_estimator=Perceptron(),
                                        n_estimators='3')
    iris = load_iris()
    assert_raise_message(ValueError,
                         "n_estimators must be an integer",
                         string_ensemble.fit, iris.data, iris.target)
    float_ensemble = BaggingClassifier(base_estimator=Perceptron(),
                                       n_estimators=3.0)
    assert_raise_message(ValueError,
                         "n_estimators must be an integer",
                         float_ensemble.fit, iris.data, iris.target)
