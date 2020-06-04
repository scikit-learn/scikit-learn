# License: BSD 3 clause

import warnings

import pytest

from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.pipeline import make_pipeline
from sklearn.compose import make_column_transformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.base import is_classifier, is_regressor, ClusterMixin, clone
from sklearn._callbacks import BaseCallback, _check_callback_params
from sklearn.utils import all_estimators
from sklearn.utils._testing import set_random_state
from sklearn.utils.estimator_checks import (
    _construct_instance,
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
    estimator = _construct_instance(Estimator)

    tags = estimator._get_tags()

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


def check_has_callback(est, callback):
    assert hasattr(est, "_callbacks") and est._callbacks is not None
    assert est._callbacks[0] is callback
    return True


def test_set_callbacks_clone():
    # Check that clone preserves callbacks
    est = StandardScaler()
    callback = CheckCallback()
    est._set_callbacks(callback)
    check_has_callback(est, callback)

    est2 = clone(est)
    check_has_callback(est2, callback)


def test_set_callbacks():
    # Check that callbacks are set recursively for meta-estimators

    X, y = load_iris(return_X_y=True)

    # check simple pipeline (recursive)
    callback = CheckCallback()
    pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=3))
    pipe._set_callbacks(callback)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        pipe.fit(X, y)
    check_has_callback(pipe, callback)
    check_has_callback(pipe.named_steps['standardscaler'], callback)
    check_has_callback(pipe.named_steps['logisticregression'], callback)

    # check simple pipeline (non recursive)
    callback = CheckCallback()
    pipe = make_pipeline(StandardScaler())
    pipe._set_callbacks(callback, deep=False)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        pipe.fit(X, y)
    check_has_callback(pipe, callback)
    assert not hasattr(pipe.named_steps['standardscaler'], "_callbacks")

    # check column transformer
    callback = CheckCallback()
    pipe = make_column_transformer(
            (StandardScaler(), [0, 1]), (MinMaxScaler(), [2, 3]),
    )

    pipe._set_callbacks(callback)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        pipe.fit(X, y)
    check_has_callback(pipe, callback)
    check_has_callback(pipe.named_transformers_['standardscaler'], callback)
    check_has_callback(pipe.named_transformers_['minmaxscaler'], callback)
