# License: BSD 3 clause

from pathlib import Path
import pytest

from sklearn.callback.tests._utils import TestingCallback
from sklearn.callback.tests._utils import TestingAutoPropagatedCallback
from sklearn.callback.tests._utils import NotValidCallback
from sklearn.callback.tests._utils import Estimator
from sklearn.callback.tests._utils import MetaEstimator


@pytest.mark.parametrize("callbacks",
    [
        TestingCallback(),
        [TestingCallback()],
        [TestingCallback(), TestingAutoPropagatedCallback()],
    ]
)
def test_set_callbacks(callbacks):
    """Sanity check for the _set_callbacks method"""
    estimator = Estimator()

    set_callbacks_return = estimator._set_callbacks(callbacks)
    assert hasattr(estimator, "_callbacks")
    assert estimator._callbacks in (callbacks, [callbacks])
    assert set_callbacks_return is estimator


@pytest.mark.parametrize("callbacks", [None, NotValidCallback()])
def test_set_callbacks_error(callbacks):
    """Check the error message when not passing a valid callback to _set_callbacks"""
    estimator = Estimator()

    with pytest.raises(TypeError, match="callbacks must be subclasses of BaseCallback"):
        estimator._set_callbacks(callbacks)


def test_propagate_callbacks():
    """Sanity check for the _propagate_callbacks method"""
    not_propagated_callback = TestingCallback()
    propagated_callback = TestingAutoPropagatedCallback()

    estimator = Estimator()
    estimator._set_callbacks([not_propagated_callback, propagated_callback])

    sub_estimator = Estimator()
    estimator._propagate_callbacks(sub_estimator, parent_node=None)

    assert hasattr(sub_estimator, "_parent_ct_node")
    assert not_propagated_callback not in sub_estimator._callbacks
    assert propagated_callback in sub_estimator._callbacks 


def test_propagate_callback_no_callback():
    """Check that no callback is propagated if there's no callback"""
    estimator = Estimator()
    sub_estimator = Estimator()
    estimator._propagate_callbacks(sub_estimator, parent_node=None)

    assert not hasattr(estimator, "_callbacks")
    assert not hasattr(sub_estimator, "_callbacks")


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
    """Check that _eval_callbacks_on_fit_begin creates and dumps the computation tree"""
    estimator = Estimator()._set_callbacks(TestingCallback())
    assert not hasattr(estimator, "_computation_tree")

    levels = [
        {"descr": "fit", "max_iter": 10},
        {"descr": "iter", "max_iter": None},
    ]
    ct_root = estimator._eval_callbacks_on_fit_begin(levels=levels)
    assert hasattr(estimator, "_computation_tree")
    assert ct_root is estimator._computation_tree.root

    ct_pickle = Path(estimator._computation_tree.tree_dir) / "computation_tree.pkl"
    assert ct_pickle.exists()


def test_callback_context_finalize():
    """Check that the folder containing the computation tree of the estimator is
    deleted when there are no reference left to its callbacks.
    """
    callback = TestingCallback()

    # estimator is not fitted, its computation tree is not built yet
    est = Estimator()._set_callbacks(callbacks=callback)
    assert not hasattr(est, "_computation_tree")

    # estimator is fitted, a folder has been created to hold its computation tree
    est.fit(X=None, y=None)
    assert hasattr(est, "_computation_tree")
    tree_dir = est._computation_tree.tree_dir
    assert tree_dir.is_dir()

    # there is no more reference to the estimator, but there is still a reference to the
    # callback which might need to access the computation tree
    del est
    assert tree_dir.is_dir()

    # there is no more reference to the callback, the computation tree folder must be
    # deleted
    del callback
    assert not tree_dir.is_dir()
