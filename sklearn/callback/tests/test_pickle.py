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
