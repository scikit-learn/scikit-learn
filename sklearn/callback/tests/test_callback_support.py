# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn.base import clone
from sklearn.callback._callback_support import init_callback_context
from sklearn.callback.tests._utils import (
    FailingCallback,
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


@pytest.mark.parametrize("fail_at", ["on_fit_begin", "on_fit_task_end", "on_fit_end"])
def test_callback_error(fail_at):
    """Check that a failing callback is properly teared down."""
    callback = FailingCallback(fail_at=fail_at)
    estimator = MaxIterEstimator().set_callbacks(callback)
    with pytest.raises(ValueError, match=f"Failing callback failed at {fail_at}"):
        estimator.fit()

    assert callback.count_hooks("on_fit_begin") == 1
    assert callback.count_hooks("on_fit_end") == 1


def _fit_one(estimator, caller, context):
    """Clone an estimator and fit it."""
    est = clone(estimator)
    context.propagate_callbacks(est)
    est.fit()

    context.eval_on_function_task_end(caller)


def my_function(estimator, n_fits=4, n_jobs=2, callbacks=None):
    """Run clone+fit on the estimator in a parallel loop."""
    caller = "_fit_estimator_in_parallel"

    try:
        context = init_callback_context(
            caller, callbacks, task_name="run", max_subtasks=n_fits
        ).eval_on_function_begin(caller)

        Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_fit_one)(
                estimator,
                caller,
                context=context.subcontext(task_name="fit estimator", task_id=i),
            )
            for i in range(n_fits)
        )
    finally:
        context.eval_on_function_end(caller)


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_callbacks_with_function(n_jobs):
    """Check that callbacks work when an estimator with a callback is passed to a
    function that runs a parallel loop where each step calls a function that clones
    and fits the estimator.
    """
    n_fits, max_iter = 4, 10
    callback = TestingCallback()
    estimator = MaxIterEstimator(
        max_iter=max_iter, computation_intensity=0
    ).set_callbacks(callback)
    my_function(estimator, n_fits=n_fits, n_jobs=n_jobs)

    assert callback.count_hooks("on_fit_begin") == n_fits
    assert callback.count_hooks("on_fit_task_end") == n_fits * max_iter
    assert callback.count_hooks("on_fit_end") == n_fits
    assert callback.count_hooks("on_function_begin") == 0
    assert callback.count_hooks("on_function_task_end") == 0
    assert callback.count_hooks("on_function_end") == 0


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_autopropagated_callbacks_with_function(n_jobs):
    """Check that callbacks work when an estimator with a callback is passed to a
    function that runs a parallel loop where each step calls a function that clones
    and fits the estimator.
    """
    n_fits, max_iter = 4, 10
    callback = TestingAutoPropagatedCallback()
    estimator = MaxIterEstimator(max_iter=max_iter, computation_intensity=0)
    my_function(estimator, n_fits=n_fits, n_jobs=n_jobs, callbacks=callback)

    assert callback.count_hooks("on_fit_begin") == 0
    assert callback.count_hooks("on_fit_task_end") == n_fits * max_iter
    assert callback.count_hooks("on_fit_end") == 0
    assert callback.count_hooks("on_function_begin") == 1
    assert callback.count_hooks("on_function_task_end") == n_fits
    assert callback.count_hooks("on_function_end") == 1
