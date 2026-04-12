# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import textwrap

import pytest

import sklearn.callback._callback_support as callback_support
from sklearn.base import clone
from sklearn.callback.tests._utils import (
    FailingCallback,
    MaxIterEstimator,
    NotValidCallback,
    TestingAutoPropagatedCallback,
    TestingCallback,
)
from sklearn.utils._testing import assert_run_python_script_without_output
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


@pytest.mark.parametrize("callbacks", [None, NotValidCallback(), TestingCallback])
def test_set_callbacks_error(callbacks):
    """Check the error message when not passing a valid callback to `set_callbacks`."""
    estimator = MaxIterEstimator()

    with pytest.raises(
        TypeError,
        match="callbacks must be instances following the FitCallback protocol.",
    ):
        estimator.set_callbacks(callbacks)


@pytest.mark.parametrize(
    "fail_at", ["setup", "on_fit_task_begin", "on_fit_task_end", "teardown"]
)
def test_callback_error(fail_at):
    """Check that a failing callback is properly teared down."""
    callback = FailingCallback(fail_at=fail_at)
    estimator = MaxIterEstimator().set_callbacks(callback)
    with pytest.raises(ValueError, match=f"Failing callback failed at {fail_at}"):
        estimator.fit()

    assert callback.count_hooks("setup") == 1
    assert callback.count_hooks("teardown") == 1


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("prefer", ["threads", "processes"])
@pytest.mark.parametrize("Callback", [TestingCallback, TestingAutoPropagatedCallback])
def test_function_no_callback_support(n_jobs, prefer, Callback):
    """Check callbacks on estimators within function not supporting callbacks.

    Since the outer function does not support callbacks, there's no shared root context
    and the context trees of each sub-estimator are independent. As a result, the
    callback acts as a regular non-propagated callback: its on_fit_begin and on_fit_end
    are called once for each fit of the sub-estimator and the number of tasks is the sum
    of the number of tasks from all the sub-estimators.
    """

    def clone_and_fit(estimator):
        clone(estimator).fit()

    def func(estimator, n_fits, n_jobs, prefer):
        Parallel(n_jobs=n_jobs, prefer=prefer)(
            delayed(clone_and_fit)(estimator) for _ in range(n_fits)
        )

    n_fits, max_iter = 5, 7
    callback = Callback()
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)

    func(estimator, n_fits, n_jobs, prefer)

    assert callback.count_hooks("setup") == n_fits
    # 1 root + max_iter leaves per fit
    assert callback.count_hooks("on_fit_task_begin") == n_fits * (1 + max_iter)
    assert callback.count_hooks("on_fit_task_end") == n_fits * (1 + max_iter)
    assert callback.count_hooks("teardown") == n_fits


def test_instantiate_manager_outside_main_module():
    """Test instantiating the callbacks manager outside the __main__ module."""
    code = """
    from sklearn.callback._callback_support import get_callback_manager

    get_callback_manager()
    """
    assert_run_python_script_without_output(textwrap.dedent(code))


def test_is_debugger_session_active_requires_trace_and_debug_module(monkeypatch):
    monkeypatch.setattr(callback_support.sys, "gettrace", lambda: object())
    monkeypatch.setitem(callback_support.sys.modules, "debugpy", object())

    assert callback_support._is_debugger_session_active()

    monkeypatch.setattr(callback_support.sys, "gettrace", lambda: None)
    assert not callback_support._is_debugger_session_active()

    monkeypatch.setattr(callback_support.sys, "gettrace", lambda: object())
    monkeypatch.delitem(callback_support.sys.modules, "debugpy", raising=False)
    monkeypatch.delitem(callback_support.sys.modules, "pydevd", raising=False)
    monkeypatch.delitem(callback_support.sys.modules, "_pydevd_bundle", raising=False)
    assert not callback_support._is_debugger_session_active()


def test_is_debugger_session_active_with_pydevd_global_debugger(monkeypatch):
    class DummyPydevd:
        @staticmethod
        def get_global_debugger():
            return object()

    monkeypatch.setattr(callback_support.sys, "gettrace", lambda: None)
    monkeypatch.setitem(callback_support.sys.modules, "pydevd", DummyPydevd())

    assert callback_support._is_debugger_session_active()


def test_get_fork_context_returns_none_on_non_posix(monkeypatch):
    monkeypatch.setattr(callback_support.os, "name", "nt")
    monkeypatch.setattr(
        callback_support,
        "get_multiprocessing_context",
        lambda *_: pytest.fail("fork context should not be requested"),
    )

    assert callback_support._get_fork_context() is None


def test_get_fork_context_returns_none_when_fork_is_unavailable(monkeypatch):
    monkeypatch.setattr(callback_support.os, "name", "posix")

    def raise_value_error(name):
        assert name == "fork"
        raise ValueError

    monkeypatch.setattr(
        callback_support, "get_multiprocessing_context", raise_value_error
    )

    assert callback_support._get_fork_context() is None


def test_get_fork_context_returns_context(monkeypatch):
    sentinel = object()

    monkeypatch.setattr(callback_support.os, "name", "posix")
    monkeypatch.setattr(
        callback_support,
        "get_multiprocessing_context",
        lambda name: sentinel,
    )

    assert callback_support._get_fork_context() is sentinel


def test_create_callback_manager_prefers_loky(monkeypatch):
    calls = []

    class DummyLokyContext:
        def Manager(self):
            calls.append("loky")
            return "loky-manager"

    monkeypatch.setattr(
        callback_support, "get_loky_context", lambda: DummyLokyContext()
    )
    monkeypatch.setattr(
        callback_support,
        "_get_fork_context",
        lambda: pytest.fail("fork fallback should not be used"),
    )

    assert callback_support._create_callback_manager() == "loky-manager"
    assert calls == ["loky"]


def test_create_callback_manager_reraises_outside_debugger(monkeypatch):
    error = EOFError("loky bootstrap failed")

    class DummyLokyContext:
        def Manager(self):
            raise error

    monkeypatch.setattr(
        callback_support, "get_loky_context", lambda: DummyLokyContext()
    )
    monkeypatch.setattr(callback_support, "_is_debugger_session_active", lambda: False)
    monkeypatch.setattr(
        callback_support,
        "_get_fork_context",
        lambda: pytest.fail("fork fallback should not be used"),
    )

    with pytest.raises(EOFError) as exc_info:
        callback_support._create_callback_manager()

    assert exc_info.value is error


def test_create_callback_manager_retries_with_fork_in_debugger(monkeypatch):
    calls = []

    class DummyLokyContext:
        def Manager(self):
            calls.append("loky")
            raise EOFError("loky bootstrap failed")

    class DummyForkContext:
        def Manager(self):
            calls.append("fork")
            return "fork-manager"

    monkeypatch.setattr(
        callback_support, "get_loky_context", lambda: DummyLokyContext()
    )
    monkeypatch.setattr(callback_support, "_is_debugger_session_active", lambda: True)
    monkeypatch.setattr(
        callback_support, "_get_fork_context", lambda: DummyForkContext()
    )

    assert callback_support._create_callback_manager() == "fork-manager"
    assert calls == ["loky", "fork"]


def test_create_callback_manager_reraises_without_fork_context(monkeypatch):
    error = EOFError("loky bootstrap failed")

    class DummyLokyContext:
        def Manager(self):
            raise error

    monkeypatch.setattr(
        callback_support, "get_loky_context", lambda: DummyLokyContext()
    )
    monkeypatch.setattr(callback_support, "_is_debugger_session_active", lambda: True)
    monkeypatch.setattr(callback_support, "_get_fork_context", lambda: None)

    with pytest.raises(EOFError) as exc_info:
        callback_support._create_callback_manager()

    assert exc_info.value is error
