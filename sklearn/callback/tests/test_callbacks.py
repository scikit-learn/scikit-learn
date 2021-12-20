# License: BSD 3 clause

import pickle
import pytest
import tempfile
from time import sleep

from joblib import Parallel, delayed

from sklearn.base import BaseEstimator, clone
from sklearn.callback import Snapshot
from sklearn.callback._base import _eval_callbacks_on_fit_iter_end
from sklearn.datasets import make_classification


class Estimator(BaseEstimator):
    def __init__(self, max_iter=20):
        self.max_iter = max_iter

    def fit(self, X, y):
        root = self._eval_callbacks_on_fit_begin(
            levels=[
                {"descr": "fit", "max_iter": self.max_iter},
                {"descr": "iter", "max_iter": None},
            ],
            X=X,
            y=y,
        )

        for i in range(self.max_iter):
            if _eval_callbacks_on_fit_iter_end(
                estimator=self,
                node=root.children[i],
                reconstruction_attributes=lambda: {"n_iter_": i + 1},
            ):
                break

        self.n_iter_ = i + 1

        self._eval_callbacks_on_fit_end()

        return self


class MetaEstimator(BaseEstimator):
    def __init__(
        self, estimator, n_outer=4, n_inner=3, n_jobs=None, prefer="processes"
    ):
        self.estimator = estimator
        self.n_outer = n_outer
        self.n_inner = n_inner
        self.n_jobs = n_jobs
        self.prefer = prefer

    def fit(self, X, y):
        root = self._eval_callbacks_on_fit_begin(
            levels=[
                {"descr": "fit", "max_iter": self.n_outer},
                {"descr": "outer", "max_iter": self.n_inner},
                {"descr": "inner", "max_iter": None},
            ],
            X=X,
            y=y,
        )

        res = Parallel(n_jobs=self.n_jobs, prefer=self.prefer)(
            delayed(self._func)(self.estimator, X, y, node, i)
            for i, node in enumerate(root.children)
        )

        self._eval_callbacks_on_fit_end()

        return self

    def _func(self, estimator, X, y, parent_node, i):
        for j, node in enumerate(parent_node.children):
            est = clone(estimator)
            self._propagate_callbacks(est, parent_node=node)
            est.fit(X, y)

            _eval_callbacks_on_fit_iter_end(estimator=self, node=node)

        _eval_callbacks_on_fit_iter_end(estimator=self, node=parent_node)

        return


@pytest.mark.parametrize("n_jobs", (1, 2))
@pytest.mark.parametrize("prefer", ("threads", "processes"))
def test_snapshot_meta_estimator(n_jobs, prefer):
    # Test for the Snapshot callback
    X, y = make_classification()
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
