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
    """Sanity check for the `_set_callbacks` method."""
    estimator = Estimator()

    set_callbacks_return = estimator._set_callbacks(callbacks)
    assert hasattr(estimator, "_skl_callbacks")

    expected_callbacks = [callbacks] if not isinstance(callbacks, list) else callbacks
    assert estimator._skl_callbacks == expected_callbacks

    assert set_callbacks_return is estimator


@pytest.mark.parametrize("callbacks", [None, NotValidCallback()])
def test_set_callbacks_error(callbacks):
    """Check the error message when not passing a valid callback to `_set_callbacks`."""
    estimator = Estimator()

    with pytest.raises(TypeError, match="callbacks must be subclasses of BaseCallback"):
        estimator._set_callbacks(callbacks)


def test_propagate_callbacks():
    """Sanity check for the `_propagate_callbacks` method."""
    not_propagated_callback = TestingCallback()
    propagated_callback = TestingAutoPropagatedCallback()

    estimator = Estimator()
    metaestimator = MetaEstimator(estimator)
    metaestimator._set_callbacks([not_propagated_callback, propagated_callback])

    metaestimator._propagate_callbacks(estimator, parent_node=None)

    assert hasattr(estimator, "_parent_node")
    assert not_propagated_callback not in estimator._skl_callbacks
    assert propagated_callback in estimator._skl_callbacks


def test_propagate_callback_no_callback():
    """Check that no callback is propagated if there's no callback."""
    estimator = Estimator()
    metaestimator = MetaEstimator(estimator)
    metaestimator._propagate_callbacks(estimator, parent_node=None)

    assert not hasattr(metaestimator, "_skl_callbacks")
    assert not hasattr(estimator, "_skl_callbacks")


def test_auto_propagated_callbacks():
    """Check that it's not possible to set an auto-propagated callback on the
    sub-estimator of a meta-estimator.
    """
    estimator = Estimator()
    estimator._set_callbacks(TestingAutoPropagatedCallback())

    meta_estimator = MetaEstimator(estimator=estimator)

    match = (
        r"sub-estimators .*of a meta-estimator .*can't have auto-propagated callbacks"
    )
    with pytest.raises(TypeError, match=match):
        meta_estimator.fit(X=None, y=None)


def test_eval_callbacks_on_fit_begin():
    """Check that `_eval_callbacks_on_fit_begin` creates the computation tree."""
    estimator = Estimator()._set_callbacks(TestingCallback())
    assert not hasattr(estimator, "_computation_tree")

    tree_structure = [
        {"stage": "fit", "n_children": 10},
        {"stage": "iter", "n_children": None},
    ]
    estimator._eval_callbacks_on_fit_begin(tree_structure=tree_structure, data={})
    assert hasattr(estimator, "_computation_tree")


def test_no_callback_early_stop():
    """Check that `eval_callbacks_on_fit_iter_end` doesn't trigger early stopping
    when there's no callback.
    """
    estimator = Estimator()
    estimator.fit(X=None, y=None)

    assert estimator.n_iter_ == estimator.max_iter
