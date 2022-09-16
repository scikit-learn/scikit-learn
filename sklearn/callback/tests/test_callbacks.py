# License: BSD 3 clause

import pickle
import pytest
import tempfile

import numpy as np

from sklearn.callback import ConvergenceMonitor
from sklearn.callback import EarlyStopping
from sklearn.callback import ProgressBar
from sklearn.callback import Snapshot
from sklearn.callback import TextVerbose
from sklearn.callback.tests._utils import Estimator
from sklearn.callback.tests._utils import MetaEstimator


X = np.zeros((100, 3))
y = np.zeros(100, dtype=int)


@pytest.mark.parametrize("Callback", [ConvergenceMonitor, EarlyStopping, ProgressBar, Snapshot, TextVerbose,])
def test_callback_doesnt_hold_ref_to_estimator(Callback):
    callback = Callback()
    est = Estimator()._set_callbacks(callbacks=callback)
    est.fit(X, y)

    tree_dir = est._computation_tree.tree_dir

    del est
    del callback
    assert not tree_dir.is_dir()


@pytest.mark.parametrize("n_jobs", (1, 2))
@pytest.mark.parametrize("prefer", ("threads", "processes"))
def test_snapshot_meta_estimator(n_jobs, prefer):
    """Test for the Snapshot callback"""
    estimator = Estimator(max_iter=20)

    with tempfile.TemporaryDirectory() as tmp_dir:
        keep_last_n = 5
        callback = Snapshot(keep_last_n=keep_last_n, base_dir=tmp_dir)
        estimator._set_callbacks(callback)
        metaestimator = MetaEstimator(
            estimator=estimator, n_outer=4, n_inner=3, n_jobs=n_jobs, prefer=prefer
        )

        metaestimator.fit(X, y)

        # There's a subdir of base_dir for each clone of estimator fitted in
        # metaestimator. There are n_outer * n_inner such clones
        snapshot_dirs = list(callback.base_dir.iterdir())
        assert len(snapshot_dirs) == metaestimator.n_outer * metaestimator.n_inner

        for snapshot_dir in snapshot_dirs:
            snapshots = sorted(snapshot_dir.iterdir())
            assert len(snapshots) == keep_last_n

            for i, snapshot in enumerate(snapshots):
                with open(snapshot, "rb") as f:
                    loaded_estimator = pickle.load(f)

                # We kept last 5 snapshots out of 20 iterations.
                # This one is the 16 + i-th.
                assert loaded_estimator.n_iter_ == 16 + i


