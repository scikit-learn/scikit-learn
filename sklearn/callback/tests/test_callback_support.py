# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

import pytest

from sklearn.base import clone
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    NotValidCallback,
    TestingAutoPropagatedCallback,
    TestingCallback,
)
from sklearn.utils.parallel import Parallel, delayed


@pytest.mark.parametrize(
    "callbacks",
    [
        TestingCallback(),
        [TestingCallback()],
        [TestingCallback(), TestingAutoPropagatedCallback()],
    ],
)
def test_set_callbacks(callbacks):
    """Sanity check for the `set_callbacks` method."""
    estimator = MaxIterEstimator()

    set_callbacks_return = estimator.set_callbacks(callbacks)
    assert hasattr(estimator, "_skl_callbacks")

    expected_callbacks = [callbacks] if not isinstance(callbacks, list) else callbacks
    assert estimator._skl_callbacks == expected_callbacks

    assert set_callbacks_return is estimator


@pytest.mark.parametrize("callbacks", [None, NotValidCallback()])
def test_set_callbacks_error(callbacks):
    """Check the error message when not passing a valid callback to `set_callbacks`."""
    estimator = MaxIterEstimator()

    with pytest.raises(TypeError, match="callbacks must follow the Callback protocol."):
        estimator.set_callbacks(callbacks)


@pytest.mark.parametrize("keep_callbacks", ["warn", True, False])
def test_clone_with_callbacks(keep_callbacks):
    """Test cloning an estimator with callbacks set on it."""
    estimator = MaxIterEstimator()
    estimator.set_callbacks(TestingCallback())

    if keep_callbacks == "warn":
        with pytest.warns(
            UserWarning, match="There are callbacks set on the estimator "
        ):
            cloned_estimator = clone(estimator)
        assert not hasattr(cloned_estimator, "_skl_callbacks")
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "error", message="There are callbacks set on the estimator "
            )
            cloned_estimator = clone(estimator, keep_callbacks=keep_callbacks)
        assert hasattr(cloned_estimator, "_skl_callbacks") == keep_callbacks


def test_function_without_callback_support():
    """Check function without callback support behavior."""
    callback = TestingCallback()
    estimator = MaxIterEstimator().set_callbacks(callback)

    def fit_one(estimator):
        estimator = clone(estimator)
        estimator.fit()

    def my_function(estimator, n_fits=4):
        """A function cloning and fitting an estimator in parallel."""
        Parallel(n_jobs=4)(delayed(fit_one)(estimator) for _ in range(n_fits))

    with pytest.warns(UserWarning, match="There are callbacks set on the estimator"):
        my_function(estimator)

    for hook in ("on_fit_begin", "on_fit_task_end", "on_fit_end"):
        assert callback.count_hooks(hook) == 0
