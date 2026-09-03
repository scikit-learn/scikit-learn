# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Common pickle round-trip tests for callbacks.

These tests guard the contract that callbacks (and estimators they are attached to)
must be picklable, and that an estimator pickled after a successful fit can be
unpickled in a fresh Python interpreter.
"""

import functools
import io
import pickle
import re
import subprocess
import sys
import textwrap

import pytest

from sklearn.callback import ProgressBar, ScoringMonitor
from sklearn.callback._transport import _message_consumers
from sklearn.callback.tests._common.estimators import MaxIterEstimator
from sklearn.datasets import make_regression


def _pb():
    pytest.importorskip("rich")
    return ProgressBar()


def _sm():
    return ScoringMonitor(scoring_train="r2")


CALLBACK_FACTORIES = [
    pytest.param(_pb, id="ProgressBar"),
    pytest.param(_sm, id="ScoringMonitor"),
]


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_estimator_with_callback_pickle_roundtrip_pre_fit(factory):
    """An estimator with the callback registered but not yet fitted is picklable."""
    estimator = MaxIterEstimator().set_callbacks(factory())
    restored = pickle.loads(pickle.dumps(estimator))
    assert type(restored) is type(estimator)
    assert len(restored._skl_callbacks) == 1


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_estimator_with_callback_pickle_roundtrip_post_fit(factory):
    """An estimator with the callback registered and fitted is picklable."""
    callback = factory()
    estimator = MaxIterEstimator(max_iter=3).set_callbacks(callback)
    estimator.fit()
    restored = pickle.loads(pickle.dumps(estimator))
    assert type(restored) is type(estimator)
    assert len(restored._skl_callbacks) == 1


def test_callbacks_refit_after_pickle_in_same_process(capsys):
    """An estimator with callbacks survives an in-process pickle round-trip.

    It also supports re-fitting after being unpickled and the callbacks accumulate new
    data from the re-fit.
    """
    pytest.importorskip("rich")

    X, y = make_regression(n_samples=30, n_features=2, random_state=0)

    sm = ScoringMonitor(scoring_train="r2")
    estimator = MaxIterEstimator(max_iter=3).set_callbacks(ProgressBar(), sm)
    estimator.fit(X=X, y=y)

    captured = capsys.readouterr()
    assert re.search(r"MaxIterEstimator - fit", captured.out)
    assert re.search(r"100%", captured.out)

    original_logs = sm.get_logs(select="all")
    assert len(original_logs) == 1

    restored = pickle.loads(pickle.dumps(estimator))
    restored.fit(X=X, y=y)

    captured = capsys.readouterr()
    assert re.search(r"MaxIterEstimator - fit", captured.out)
    assert re.search(r"100%", captured.out)

    restored_logs = restored._skl_callbacks[1].get_logs(select="all")
    assert len(restored_logs) == 2
    assert restored_logs[0].train_scores == original_logs[0].train_scores


def test_callbacks_refit_after_load_in_fresh_process(tmp_path, capsys):
    """An estimator with callbacks survives unpickling in a fresh interpreter.

    It also supports re-fitting after being unpickled and the callbacks accumulate new
    data from the re-fit.
    """
    pytest.importorskip("rich")

    X, y = make_regression(n_samples=20, n_features=3, random_state=0)

    sm = ScoringMonitor(scoring_train="r2")
    estimator = MaxIterEstimator(max_iter=3).set_callbacks(ProgressBar(), sm)
    estimator.fit(X=X, y=y)

    captured = capsys.readouterr()
    assert re.search(r"MaxIterEstimator - fit", captured.out)
    assert re.search(r"100%", captured.out)

    original_logs = sm.get_logs(select="all")
    assert len(original_logs) == 1

    pkl_path = tmp_path / "est.pkl"
    with open(pkl_path, "wb") as f:
        pickle.dump(estimator, f)

    load_script = textwrap.dedent(
        f"""
        import pickle
        from sklearn.callback import ScoringMonitor
        from sklearn.datasets import make_regression

        with open({str(pkl_path)!r}, "rb") as f:
            est = pickle.load(f)

        X, y = make_regression(n_samples=20, n_features=3, random_state=1)
        est.fit(X=X, y=y)

        restored_logs = est._skl_callbacks[1].get_logs(select="all")
        assert len(restored_logs) == 2
        assert restored_logs[0].train_scores == {original_logs[0].train_scores}
        """
    )

    result = subprocess.run(
        [sys.executable, "-c", load_script], capture_output=True, timeout=120
    )

    stdout = result.stdout.decode()
    assert re.search(r"MaxIterEstimator - fit", stdout)
    assert re.search(r"100%", stdout)


class _TraversalRecorder(pickle.Pickler):
    """A pickler that records the id of every object it walks through.

    `persistent_id` is called for every object the pickler encounters, so an object is
    recorded whichever path leads to it, not only when it is a direct attribute of the
    object being pickled,
    see https://docs.python.org/3/library/pickle.html#pickle.Pickler.persistent_id.
    """

    def __init__(self):
        super().__init__(io.BytesIO(), protocol=pickle.HIGHEST_PROTOCOL)
        self.walked_through = set()

    def persistent_id(self, obj):
        self.walked_through.add(id(obj))
        return None  # pickle `obj` as usual


def _checked(hook):
    """Wrap a callback hook so that it first checks the callback it is called on.

    The check is that the pickler does not walk through the objects that listener
    threads mutate. They are found through the registered consumers: a consumer is
    normally a method bound to the container it fills, e.g. `self._log.append` for
    ScoringMonitor or `queue.put` for ProgressBar, so the container is what it is bound
    to. A consumer bound to nothing, e.g. a closure, is skipped, since there is then no
    way to tell what it mutates.
    """

    # preserve the signature of the hook because callbacks are validated against it
    @functools.wraps(hook)
    def checked_hook(self, *args, **kwargs):
        # snapshot because another thread may register a listener concurrently
        consumers = list(_message_consumers.values())
        # the containers are held, not just their ids, which could be reused once freed
        watched = [c.__self__ for c in consumers if hasattr(c, "__self__")]

        recorder = _TraversalRecorder()
        recorder.dump(self)

        offenders = recorder.walked_through & {id(container) for container in watched}
        assert not offenders, (
            f"Pickling {self.__class__.__name__} walks through a container that a"
            " listener thread mutates concurrently. Keep it away from the pickler,"
            " either by handing over a copy of it in the callback's `__getstate__`, or"
            " by storing it outside of the callback instance."
        )

        return hook(self, *args, **kwargs)

    return checked_hook


@pytest.mark.parametrize("factory", CALLBACK_FACTORIES)
def test_listener_state_is_not_walked_by_the_pickler(factory, monkeypatch):
    """Check that pickling a callback never traverses state its listener mutates.

    An estimator carrying a callback can be pickled by a background thread, e.g. loky's
    queue feeder dispatching a task to a worker, while the listener thread of that same
    callback mutates the callback's state as messages come in. Pickling a container that
    another thread mutates breaks the dump, which joblib reports as "Could not pickle
    the task to send it to the workers".

    The check runs from a hook, i.e. while the listeners are up, which is when such a
    dispatch would happen.
    """
    callback = factory()

    for hook_name in ("on_fit_task_begin", "on_fit_task_end"):
        hook = getattr(callback.__class__, hook_name)
        monkeypatch.setattr(callback.__class__, hook_name, _checked(hook))

    MaxIterEstimator(max_iter=3).set_callbacks(callback).fit()
