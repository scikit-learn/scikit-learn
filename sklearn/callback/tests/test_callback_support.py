# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn.base import clone
from sklearn.callback.tests._utils import (
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


def test_clone_warning():
    """Test the warning when cloning an estimator with callback registered to it."""
    estimator = MaxIterEstimator()
    estimator.set_callbacks(TestingCallback())
    with pytest.warns(UserWarning, match="There are callbacks set on the estimator "):
        clone(estimator)
