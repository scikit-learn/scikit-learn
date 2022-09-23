# License: BSD 3 clause

from datetime import datetime
from pathlib import Path
import pickle

from . import BaseCallback


class Snapshot(BaseCallback):
    """Take regular snapshots of an estimator

    Parameters
    ----------
    keep_last_n : int or None, default=1
        Only the last `keep_last_n` snapshots are kept on the disk. None means all
        snapshots are kept.

    base_dir : str or pathlib.Path instance, default=None
        The directory where the snapshots should be stored. If None, they are stored in
        the current directory.
    """

    request_from_reconstruction_attributes = True

    def __init__(self, keep_last_n=1, base_dir=None):
        self.keep_last_n = keep_last_n
        if keep_last_n is not None and keep_last_n <= 0:
            raise ValueError(
                "keep_last_n must be a positive integer, got"
                f" {self.keep_last_n} instead."
            )

        self.base_dir = Path("." if base_dir is None else base_dir)

    def on_fit_begin(self, estimator, X=None, y=None):
        subdir = self._get_subdir(estimator._computation_tree)
        subdir.mkdir()

    def on_fit_iter_end(self, *, estimator, node, **kwargs):
        new_estimator = kwargs.get("from_reconstruction_attributes", None)
        if new_estimator is None:
            return

        subdir = self._get_subdir(node.computation_tree)
        snapshot_filename = f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S-%f')}.pkl"

        with open(subdir / snapshot_filename, "wb") as f:
            pickle.dump(new_estimator, f)

        if self.keep_last_n is not None:
            for snapshot in sorted(subdir.iterdir())[: -self.keep_last_n]:
                snapshot.unlink(missing_ok=True)

    def on_fit_end(self):
        pass

    def _get_subdir(self, computation_tree):
        """Return the sub directory containing the snapshots of the estimator"""
        subdir = (
            self.base_dir
            / f"snapshots_{computation_tree.estimator_name}_{str(computation_tree.uid)}"
        )

        return subdir
