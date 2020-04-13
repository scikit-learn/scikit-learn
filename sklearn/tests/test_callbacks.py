# License: BSD 3 clause

import numpy as np
import pytest

from sklearn.base import is_classifier, is_regressor, ClusterMixin
from sklearn.datasets import load_iris
from sklearn._callbacks import BaseCallback, _check_callback_params
from sklearn.utils.testing import all_estimators, set_random_state
from sklearn.utils.estimator_checks import (
    _construct_instance,
    _safe_tags,
    _enforce_estimator_tags_y,
)


@pytest.fixture(scope="module")
def iris():
    X, y = load_iris(return_X_y=True)
    return X, y.reshape(-1, 1)


def test_check_callback_params():
    _check_callback_params(n_iter=2, max_iter=10, loss=0.1)

    msg = "Invalid callback parameters: a, must be one of n_iter.*"
    with pytest.raises(ValueError, match=msg):
        _check_callback_params(a=0)

    msg = "Invalid callback parameters: max_iter=1.0 is not of type .*int"
    with pytest.raises(ValueError, match=msg):
        _check_callback_params(max_iter=1.)


def _supported_estimators():
    for (name, Estimator) in all_estimators():
        if name.startswith("_"):
            continue
        if name in [
            # need to make appropriate 1D test data
            "IsotonicRegression"
        ]:
            continue

        if (
            is_classifier(Estimator)
            or is_regressor(Estimator)
            or issubclass(Estimator, ClusterMixin)
        ):
            yield name, Estimator


class CheckCallback(BaseCallback):
    def __init__(self):
        self.n_calls = 0
        self.n_fit_calls = 0

    def fit(self, estimator, X, y):
        self.n_fit_calls += 1

    def __call__(self, **kwargs):
        self.n_calls += 1
        _check_callback_params(**kwargs)


@pytest.mark.parametrize("name, Estimator", _supported_estimators())
def test_callback(name, Estimator, iris):

    # if name in [
    #    "DictionaryLearning",
    #    "MiniBatchDictionaryLearning",
    #    "TheilSenRegressor",
    # ]:
    #    # TODO: make these tests faster
    #    pytest.skip("Slow tests, skipping for now")

    estimator = _construct_instance(Estimator)

    tags = _safe_tags(estimator)

    callback = CheckCallback()
    estimator._set_callbacks([callback])
    if tags.get("X_types", []) == ["string"]:
        X = ["some document", "another document"]
        y = None
    else:
        X, y = iris
        y = _enforce_estimator_tags_y(estimator, y)
    set_random_state(estimator, 0)

    assert callback.n_calls == 0
    estimator.fit(X, y)
    if callback.n_fit_calls == 0:
        pytest.skip("callbacks not implemented")
