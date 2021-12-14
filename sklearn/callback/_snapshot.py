# License: BSD 3 clause

from copy import copy
from datetime import datetime
from pathlib import Path
import pickle

import numpy as np

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

    Attributes
    ----------
    directory : pathlib.Path instance
        The directory where the snapshots are saved. It's a sub-directory of `base_dir`.
    """

    request_reconstruction_attributes = True

    def __init__(self, keep_last_n=1, base_dir=None):
        self.keep_last_n = keep_last_n
        if keep_last_n is not None and keep_last_n <= 0:
            raise ValueError(
                "keep_last_n must be a positive integer, got"
                f" {self.keep_last_n} instead."
            )

        self.base_dir = Path("." if base_dir is None else base_dir)

    def on_fit_begin(self, estimator, X=None, y=None):
        self.estimator = estimator

        # Use a hash in the name of this directory to avoid name collision if several
        # clones of this estimator are fitted in parallel in a meta-estimator for
        # instance.
        dir_name = (
            "snapshots_"
            f"{self.estimator.__class__.__name__}_"
            f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S-%f')}_"
            f"{hash(self.estimator._computation_tree)}"
        )

        self.directory = self.base_dir / dir_name
        self.directory.mkdir()

    def on_fit_iter_end(self, *, node, **kwargs):
        reconstruction_attributes = kwargs.get("reconstruction_attributes", None)
        if reconstruction_attributes is None:
            return

        new_estimator = copy(self.estimator)
        for key, val in reconstruction_attributes.items():
            setattr(new_estimator, key, val)

        file_name = f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S-%f')}.pkl"
        file_path = self.directory / file_name

        with open(file_path, "wb") as f:
            pickle.dump(new_estimator, f)

        if self.keep_last_n is not None:
            for snapshot in sorted(self.directory.iterdir())[: -self.keep_last_n]:
                snapshot.unlink(missing_ok=True)

    def on_fit_end(self):
        if self.keep_last_n is not None:
            for snapshot in sorted(self.directory.iterdir())[: -self.keep_last_n]:
                snapshot.unlink()
