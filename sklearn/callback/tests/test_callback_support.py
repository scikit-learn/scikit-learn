# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

import pytest

from sklearn.base import clone
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    MetaEstimator,
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


def test_no_clone_warning_with_meta_est():
    """Test that the warning is not raised with callback-compatible meta-estimator.

    When the cloning is done inside of a callback-compatible meta-estimator, the warning
    about callbacks not being cloned should not be raised.
    """
    estimator = MaxIterEstimator()
    estimator.set_callbacks(TestingCallback())
    meta_est = MetaEstimator(estimator=estimator)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        meta_est.fit()
