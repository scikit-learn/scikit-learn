# License: BSD 3 clause
# Authors: the scikit-learn developers

import pytest

from sklearn.callback.tests._utils import (
    Estimator,
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
    estimator = Estimator()

    set_callbacks_return = estimator.set_callbacks(callbacks)
    assert hasattr(estimator, "_skl_callbacks")

    expected_callbacks = [callbacks] if not isinstance(callbacks, list) else callbacks
    assert estimator._skl_callbacks == expected_callbacks

    assert set_callbacks_return is estimator


@pytest.mark.parametrize("callbacks", [None, NotValidCallback()])
def test_set_callbacks_error(callbacks):
    """Check the error message when not passing a valid callback to `set_callbacks`."""
    estimator = Estimator()

    with pytest.raises(
        TypeError, match="callbacks must follow the CallbackProtocol protocol."
    ):
        estimator.set_callbacks(callbacks)


def test_init_callback_context():
    """Sanity check for the `init_callback_context` method."""
    estimator = Estimator()
    callback_ctx = estimator.init_callback_context()

    assert hasattr(estimator, "_callback_fit_ctx")
    assert hasattr(callback_ctx, "_callbacks")


def test_propagate_callbacks():
    """Sanity check for the `propagate_callbacks` method."""
    not_propagated_callback = TestingCallback()
    propagated_callback = TestingAutoPropagatedCallback()

    estimator = Estimator()
    metaestimator = MetaEstimator(estimator)
    metaestimator.set_callbacks([not_propagated_callback, propagated_callback])

    callback_ctx = metaestimator.init_callback_context()
    callback_ctx.propagate_callbacks(estimator)

    assert hasattr(estimator, "_parent_callback_ctx")
    assert not_propagated_callback not in estimator._skl_callbacks
    assert propagated_callback in estimator._skl_callbacks


def test_propagate_callback_no_callback():
    """Check that no callback is propagated if there's no callback."""
    estimator = Estimator()
    metaestimator = MetaEstimator(estimator)

    callback_ctx = metaestimator.init_callback_context()
    assert len(callback_ctx._callbacks) == 0

    callback_ctx.propagate_callbacks(estimator)

    assert not hasattr(metaestimator, "_skl_callbacks")
    assert not hasattr(estimator, "_skl_callbacks")


def test_auto_propagated_callbacks():
    """Check that it's not possible to set an auto-propagated callback on the
    sub-estimator of a meta-estimator.
    """
    estimator = Estimator()
    estimator.set_callbacks(TestingAutoPropagatedCallback())
    meta_estimator = MetaEstimator(estimator=estimator)

    match = (
        r"sub-estimator .*of a meta-estimator .*can't have auto-propagated callbacks"
    )
    with pytest.raises(TypeError, match=match):
        meta_estimator.fit(X=None, y=None)
