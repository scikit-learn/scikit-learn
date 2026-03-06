# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn.callback.tests._utils import (
    FailingCallback,
    MaxIterEstimator,
    NotValidCallback,
    TestingAutoPropagatedCallback,
    TestingCallback,
)


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


@pytest.mark.parametrize("fail_at", ["on_fit_begin", "on_fit_task_end", "on_fit_end"])
def test_callback_error(fail_at):
    """Check that a failing callback is properly teared down."""
    callback = FailingCallback(fail_at=fail_at)
    estimator = MaxIterEstimator().set_callbacks(callback)
    with pytest.raises(ValueError, match=f"Failing callback failed at {fail_at}"):
        estimator.fit()

    assert callback.count_hooks("on_fit_begin") == 1
    assert callback.count_hooks("on_fit_end") == 1
