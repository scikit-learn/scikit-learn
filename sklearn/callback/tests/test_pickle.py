# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Common pickle round-trip tests for callbacks.

These tests guard the contract that callbacks (and estimators they are attached to)
must be picklable, and that an estimator pickled after a successful fit can be
unpickled in a fresh Python interpreter.
"""

import pickle
import re
import subprocess
import sys
import textwrap

import pytest

from sklearn.callback import ProgressBar, ScoringMonitor
from sklearn.callback.tests._utils import MaxIterEstimator
from sklearn.datasets import make_regression


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

    sm = ScoringMonitor(scoring="r2")
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
    assert restored_logs[0].data == original_logs[0].data


def test_callbacks_refit_after_load_in_fresh_process(tmp_path, capsys):
    """An estimator with callbacks survives unpickling in a fresh interpreter.

    It also supports re-fitting after being unpickled and the callbacks accumulate new
    data from the re-fit.
    """
    pytest.importorskip("rich")

    X, y = make_regression(n_samples=20, n_features=3, random_state=0)

    sm = ScoringMonitor(scoring="r2")
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
        assert restored_logs[0].data == {original_logs[0].data}
        """
    )

    result = subprocess.run(
        [sys.executable, "-c", load_script], capture_output=True, timeout=120
    )

    stdout = result.stdout.decode()
    assert re.search(r"MaxIterEstimator - fit", stdout)
    assert re.search(r"100%", stdout)
