"""
Testing for the base module (sklearn.ensemble.base).
"""

# Authors: Gilles Louppe
# License: BSD 3 clause

import numpy as np
from numpy.testing import assert_equal

from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_true

from sklearn.datasets import load_iris
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble.base import _set_random_states
from sklearn.linear_model import Perceptron
from collections import OrderedDict
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectFromModel


def test_base():
    # Check BaseEnsemble methods.
    ensemble = BaggingClassifier(base_estimator=Perceptron(random_state=None),
                                 n_estimators=3)

    iris = load_iris()
    ensemble.fit(iris.data, iris.target)
    ensemble.estimators_ = []  # empty the list and create estimators manually

    ensemble._make_estimator()
    random_state = np.random.RandomState(3)
    ensemble._make_estimator(random_state=random_state)
    ensemble._make_estimator(random_state=random_state)
    ensemble._make_estimator(append=False)

    assert_equal(3, len(ensemble))
    assert_equal(3, len(ensemble.estimators_))

    assert_true(isinstance(ensemble[0], Perceptron))
    assert_equal(ensemble[0].random_state, None)
    assert_true(isinstance(ensemble[1].random_state, int))
    assert_true(isinstance(ensemble[2].random_state, int))
    assert_not_equal(ensemble[1].random_state, ensemble[2].random_state)

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


def test_set_random_states():
    # Linear Discriminant Analysis doesn't have random state: smoke test
    _set_random_states(LinearDiscriminantAnalysis(), random_state=17)

    clf1 = Perceptron(random_state=None)
    assert_equal(clf1.random_state, None)
    # check random_state is None still sets
    _set_random_states(clf1, None)
    assert_true(isinstance(clf1.random_state, int))

    # check random_state fixes results in consistent initialisation
    _set_random_states(clf1, 3)
    assert_true(isinstance(clf1.random_state, int))
    clf2 = Perceptron(random_state=None)
    _set_random_states(clf2, 3)
    assert_equal(clf1.random_state, clf2.random_state)

    # nested random_state

    def make_steps():
        return [('sel', SelectFromModel(Perceptron(random_state=None))),
                ('clf', Perceptron(random_state=None))]

    est1 = Pipeline(make_steps())
    _set_random_states(est1, 3)
    assert_true(isinstance(est1.steps[0][1].estimator.random_state, int))
    assert_true(isinstance(est1.steps[1][1].random_state, int))
    assert_not_equal(est1.get_params()['sel__estimator__random_state'],
                     est1.get_params()['clf__random_state'])

    # ensure multiple random_state paramaters are invariant to get_params()
    # iteration order

    class AlphaParamPipeline(Pipeline):
        def get_params(self, *args, **kwargs):
            params = Pipeline.get_params(self, *args, **kwargs).items()
            return OrderedDict(sorted(params))

    class RevParamPipeline(Pipeline):
        def get_params(self, *args, **kwargs):
            params = Pipeline.get_params(self, *args, **kwargs).items()
            return OrderedDict(sorted(params, reverse=True))

    for cls in [AlphaParamPipeline, RevParamPipeline]:
        est2 = cls(make_steps())
        _set_random_states(est2, 3)
        assert_equal(est1.get_params()['sel__estimator__random_state'],
                     est2.get_params()['sel__estimator__random_state'])
        assert_equal(est1.get_params()['clf__random_state'],
                     est2.get_params()['clf__random_state'])
