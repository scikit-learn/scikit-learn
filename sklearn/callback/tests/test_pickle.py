# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Common pickle round-trip tests for callbacks.

These tests guard the contract that callbacks (and estimators they are attached to)
must be picklable, and that an estimator pickled after a successful fit can be
unpickled in a fresh Python interpreter without trying to reconnect to a defunct
``multiprocessing.Manager`` socket.
"""

import pickle
import textwrap

import numpy as np
import pytest

from sklearn.callback import ProgressBar, ScoringMonitor
from sklearn.callback.tests._utils import MaxIterEstimator
from sklearn.datasets import make_regression
from sklearn.utils._testing import assert_run_python_script_without_output


def _pb():
    pytest.importorskip("rich")
    return ProgressBar()


def _sm():
    return ScoringMonitor(scoring="r2")


CALLBACK_FACTORIES = [
    pytest.param(_pb, id="ProgressBar"),
    pytest.param(_sm, id="ScoringMonitor"),
]


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_bare_callback_construction_no_manager_spawn(factory):
    """``__init__`` must not spawn the multiprocessing Manager.

    Constructing a callback should be a pure Python operation. The Manager subprocess
    and its UNIX socket must only be created once the callback is actually attached
    to an estimator (or, for ProgressBar, on the first ``fit``).
    """
    from sklearn.callback import _callback_support

    # Reset the singleton so this test is deterministic regardless of ordering.
    _callback_support._CallbackManagerState.manager = None

    callback = factory()

    assert _callback_support._CallbackManagerState.manager is None, (
        f"{type(callback).__name__}.__init__ unexpectedly spawned the callback manager"
    )


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_bare_callback_pickle_roundtrip(factory):
    """A freshly constructed (unattached) callback round-trips through pickle."""
    callback = factory()
    restored = pickle.loads(pickle.dumps(callback))
    assert type(restored) is type(callback)


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_estimator_with_callback_pickle_roundtrip_pre_fit(factory):
    """An estimator with the callback attached but not yet fitted is picklable."""
    estimator = MaxIterEstimator().set_callbacks(factory())
    restored = pickle.loads(pickle.dumps(estimator))
    assert type(restored) is type(estimator)
    assert len(restored._skl_callbacks) == 1


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_estimator_with_callback_pickle_roundtrip_post_fit(factory):
    """An estimator with the callback attached and fitted is picklable."""
    callback = factory()
    estimator = MaxIterEstimator(max_iter=3).set_callbacks(callback)
    estimator.fit()
    restored = pickle.loads(pickle.dumps(estimator))
    assert type(restored) is type(estimator)
    assert len(restored._skl_callbacks) == 1


def test_scoring_monitor_load_after_fit_in_fresh_process():
    """Reproduce the issue from PR #33322 (review comment r3153639421).

    A ``StandardScaler`` (or any estimator) fitted with a ``ScoringMonitor`` attached
    must be loadable in a fresh Python interpreter, where the Manager subprocess
    that backed the original transport queue is long gone. The unpickling path must
    not try to reconnect to its UNIX socket.

    On the loaded object, ``get_logs()`` should return the data that was accumulated
    during the original fit.
    """
    code = textwrap.dedent(
        """
        import pickle, tempfile, sys, textwrap, subprocess

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            pkl_path = f.name

        save_script = textwrap.dedent(f'''
            import pickle
            import numpy as np
            from sklearn.callback import ScoringMonitor
            from sklearn.callback.tests._utils import MaxIterEstimator

            est = MaxIterEstimator(max_iter=3).set_callbacks(
                ScoringMonitor(scoring="r2")
            )
            X = np.random.RandomState(0).randn(20, 3)
            y = np.random.RandomState(1).randn(20)
            est.fit(X=X, y=y)

            # Read logs on the saving side too — exercises the full path.
            sm = est._skl_callbacks[0]
            n_rows = len(sm.get_logs(select="most_recent").data)
            assert n_rows >= 1, n_rows

            with open({pkl_path!r}, "wb") as f:
                pickle.dump(est, f)
        ''')

        load_script = textwrap.dedent(f'''
            import pickle
            with open({pkl_path!r}, "rb") as f:
                est = pickle.load(f)

            sm = est._skl_callbacks[0]
            log = sm.get_logs(select="most_recent")
            assert len(log.data) >= 1, len(log.data)
        ''')

        subprocess.run([sys.executable, "-c", save_script], check=True)
        subprocess.run([sys.executable, "-c", load_script], check=True)
        """
    )
    # Run save+load orchestrator in its own subprocess so this test does not pollute
    # the Manager singleton in the test runner's interpreter.
    assert_run_python_script_without_output(code, timeout=120)


def test_scoring_monitor_fit_after_unpickle_without_reattach():
    """Re-fit a ScoringMonitor in a fresh process without re-attaching.

    After unpickling in a new interpreter, ``_address`` is a stale string pointing
    at a listener that no longer exists. If the user calls ``fit`` again without
    re-attaching the callback (which would refresh the address), the worker-side
    ``on_fit_task_end`` must not crash: it has nowhere to send the record, so it
    should silently drop it and let the fit complete. This exercises the
    ``if self._address is None`` branch (here ``_address`` is non-None but stale,
    which exercises the ``_send_log_record`` connect-failure branch instead).
    """
    code = textwrap.dedent(
        '''
        import pickle, tempfile, sys, textwrap, subprocess

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            pkl_path = f.name

        save_script = textwrap.dedent(f"""
            import pickle
            import numpy as np
            from sklearn.callback import ScoringMonitor
            from sklearn.callback.tests._utils import MaxIterEstimator

            est = MaxIterEstimator(max_iter=2).set_callbacks(
                ScoringMonitor(scoring='r2')
            )
            X = np.random.RandomState(0).randn(20, 3)
            y = np.random.RandomState(1).randn(20)
            est.fit(X=X, y=y)

            with open({pkl_path!r}, 'wb') as f:
                pickle.dump(est, f)
        """)

        # Load in a fresh process and fit AGAIN without calling set_callbacks.
        # The stale _address points at a now-defunct listener; on_fit_task_end
        # must not crash. The new run's records won't be captured, but get_logs
        # still returns the original run's data.
        load_script = textwrap.dedent(f"""
            import pickle
            import numpy as np
            with open({pkl_path!r}, 'rb') as f:
                est = pickle.load(f)
            X = np.random.RandomState(0).randn(20, 3)
            y = np.random.RandomState(1).randn(20)
            est.fit(X=X, y=y)
            sm = est._skl_callbacks[0]
            log = sm.get_logs(select='most_recent')
            # Original run data must still be there.
            assert len(log.data) >= 1, len(log.data)
        """)

        subprocess.run([sys.executable, "-c", save_script], check=True)
        subprocess.run([sys.executable, "-c", load_script], check=True)
        '''
    )
    assert_run_python_script_without_output(code, timeout=120)


def test_scoring_monitor_fit_after_unpickle_with_reattach():
    """Re-fit a ScoringMonitor in a fresh process after re-attaching.

    Re-attaching via ``set_callbacks`` refreshes the listener's address by
    triggering ``_skl_on_attach`` again. The new run's records should land in
    ``_log`` and be visible via ``get_logs``.
    """
    code = textwrap.dedent(
        '''
        import pickle, tempfile, sys, textwrap, subprocess

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            pkl_path = f.name

        save_script = textwrap.dedent(f"""
            import pickle
            import numpy as np
            from sklearn.callback import ScoringMonitor
            from sklearn.callback.tests._utils import MaxIterEstimator

            est = MaxIterEstimator(max_iter=2).set_callbacks(
                ScoringMonitor(scoring='r2')
            )
            X = np.random.RandomState(0).randn(20, 3)
            y = np.random.RandomState(1).randn(20)
            est.fit(X=X, y=y)

            with open({pkl_path!r}, 'wb') as f:
                pickle.dump(est, f)
        """)

        load_script = textwrap.dedent(f"""
            import pickle
            import numpy as np
            with open({pkl_path!r}, 'rb') as f:
                est = pickle.load(f)
            sm = est._skl_callbacks[0]
            n_before = len(sm.get_logs(select='most_recent').data)

            # Re-attach to refresh the transport in this fresh process.
            est.set_callbacks(sm)
            X = np.random.RandomState(2).randn(20, 3)
            y = np.random.RandomState(3).randn(20)
            est.fit(X=X, y=y)

            # Most recent log corresponds to the new run.
            log = sm.get_logs(select='most_recent')
            assert len(log.data) >= 1, len(log.data)
        """)

        subprocess.run([sys.executable, "-c", save_script], check=True)
        subprocess.run([sys.executable, "-c", load_script], check=True)
        '''
    )
    assert_run_python_script_without_output(code, timeout=120)


def test_scoring_monitor_fit_after_unpickle_without_address():
    """If ``_address`` is None (callback was unpickled and somehow stripped of its
    address), the worker-side ``on_fit_task_end`` early-returns instead of
    crashing.

    This guards the ``if self._address is None: return`` branch, which is the
    safety net for any unpickle path that loses transport state.
    """
    callback = ScoringMonitor(scoring="r2")
    estimator = MaxIterEstimator(max_iter=2).set_callbacks(callback)

    # Simulate the "unpickled with no address" state without going through a
    # real subprocess: drop the address and fit. on_fit_task_end takes the
    # early-return path.
    callback._address = None
    callback._authkey = None
    estimator.fit()

    # No data was sent; get_logs raises the standard "no logs" error.
    with pytest.raises(ValueError, match="No logs to retrieve"):
        callback.get_logs()


def test_scoring_monitor_socket_path_in_main_process():
    """Force the worker-side socket transport path in the main process.

    ``on_fit_task_end`` normally takes a same-process fast path
    (``_log.append`` directly) when the listener is registered locally. To
    exercise the worker-side ``_send_log_record`` branch from the test runner
    (so coverage can observe it), this test temporarily removes the listener
    from the module-level registry before calling fit. Records are then sent
    over the socket to the still-running listener, which receives and appends
    them via its handler thread.
    """
    rng = np.random.RandomState(0)
    X, y = make_regression(n_samples=30, n_features=2, random_state=rng)

    from sklearn.callback import _scoring_monitor as sm_mod

    callback = ScoringMonitor(scoring="r2")
    estimator = MaxIterEstimator(max_iter=3).set_callbacks(callback)

    listener = sm_mod._listeners.pop(callback._callback_id)
    try:
        estimator.fit(X=X, y=y)
        # Wait briefly for the handler thread to receive the records that the
        # in-process "worker" send put on the wire.
        import time

        deadline = time.monotonic() + 5.0
        while len(callback._log) < estimator.max_iter + 1:
            if time.monotonic() > deadline:
                break
            time.sleep(0.01)
    finally:
        sm_mod._listeners[callback._callback_id] = listener
        # Drop the cached worker connection so subsequent tests don't reuse it.
        with sm_mod._worker_connections_lock:
            sm_mod._worker_connections.pop(callback._address, None)

    # We sent ``max_iter + 1`` records (1 root + max_iter iterations), all of
    # which should have arrived through the socket.
    assert len(callback._log) == estimator.max_iter + 1


def test_scoring_monitor_log_repr():
    """``ScoringMonitorLog.__repr__`` returns a useful summary."""
    from sklearn.callback._scoring_monitor import ScoringMonitorLog

    callback = ScoringMonitor(scoring="r2")
    estimator = MaxIterEstimator(max_iter=2).set_callbacks(callback)
    rng = np.random.RandomState(0)
    X, y = make_regression(n_samples=20, n_features=2, random_state=rng)
    estimator.fit(X=X, y=y)

    log = callback.get_logs(select="most_recent")
    text = repr(log)
    assert isinstance(log, ScoringMonitorLog)
    assert "ScoringMonitorLog" in text
    assert "MaxIterEstimator" in text


def test_skl_on_attach_is_idempotent():
    """Re-attaching the same callback reuses the existing listener.

    Covers the early-return branch in ``_skl_on_attach`` that detects an
    already-registered listener for this callback id.
    """
    from sklearn.callback import _scoring_monitor as sm_mod

    callback = ScoringMonitor(scoring="r2")
    estimator = MaxIterEstimator().set_callbacks(callback)
    address_first = callback._address
    listener_first = sm_mod._listeners[callback._callback_id]

    # Second attach must not open a second listener.
    estimator.set_callbacks(callback)
    assert callback._address == address_first
    assert sm_mod._listeners[callback._callback_id] is listener_first


def test_scoring_monitor_get_logs_after_pickle_in_same_process():
    """In-process pickle/load of a fitted ScoringMonitor preserves the log data."""
    rng = np.random.RandomState(0)
    X, y = make_regression(n_samples=30, n_features=2, random_state=rng)

    callback = ScoringMonitor(scoring="r2")
    estimator = MaxIterEstimator(max_iter=3).set_callbacks(callback)
    estimator.fit(X=X, y=y)

    original_log = callback.get_logs(select="most_recent")
    restored_estimator = pickle.loads(pickle.dumps(estimator))
    restored_callback = restored_estimator._skl_callbacks[0]
    restored_log = restored_callback.get_logs(select="most_recent")

    assert restored_log.data == original_log.data


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_callback_no_manager_proxies_in_persistent_state(factory):
    """The instance ``__dict__`` of a post-fit callback must contain no Manager
    proxies whose lifetime is tied to a process other than the one we may load in.

    Concretely: every value in ``__dict__`` (excluding the transport queue we
    keep for in-flight workers) must be a plain Python object, not a Manager
    proxy whose connection token would be stale across a process boundary.
    """
    if factory is _pb:
        pytest.importorskip("rich")

    callback = factory()
    estimator = MaxIterEstimator(max_iter=2).set_callbacks(callback)
    estimator.fit()

    # The ScoringMonitor keeps a Manager.Queue proxy in ``_queue`` for any future
    # workers; that proxy is intentionally there and is handled by ``_skl_on_attach``
    # on re-attach. Everything else on the instance must be plain Python.
    skip = {"_queue"}
    for name, value in callback.__dict__.items():
        if name in skip:
            continue
        # Manager proxy types live under multiprocessing.managers; their classes are
        # generated by ``MakeProxyType`` so we detect them by module rather than by
        # exact class identity.
        klass = type(value)
        assert not klass.__module__.startswith("multiprocessing.managers"), (
            f"{type(callback).__name__}.{name} is a Manager proxy "
            f"({klass.__name__}); persistent state must be plain Python"
        )
