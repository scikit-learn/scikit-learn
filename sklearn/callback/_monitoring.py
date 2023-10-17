# License: BSD 3 clause

# import os
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import pandas as pd

from . import BaseCallback


class Monitoring(BaseCallback):
    """Monitor model convergence.

    Parameters
    ----------
    monitor :

    X_val : ndarray, default=None
        Validation data

    y_val : ndarray, default=None
        Validation target

    Attributes
    ----------
    data : pandas.DataFrame
        The monitored quantities at each iteration.
    """

    request_from_reconstruction_attributes = True

    def __init__(
        self,
        *,
        monitor="objective_function",
        on="validation_set",
        validation_split="auto",
    ):
        from ..model_selection import KFold

        self.validation_split = validation_split
        if validation_split == "auto":
            self.validation_split = KFold(n_splits=5, shuffle=True, random_state=42)
        self.monitor = monitor
        self.on = on

        self._data_dir = TemporaryDirectory()
        self._data_files = {}

        if isinstance(self.monitor, str):
            self.monitor_name = self.monitor
        elif callable(self.monitor):
            self.monitor_name = self.monitor.__name__

    def on_fit_begin(self, estimator, *, X=None, y=None):
        fname = Path(self._data_dir.name) / f"{estimator._computation_tree.uid}.csv"
        with open(fname, "w") as file:
            file.write(f"iteration,{self.monitor_name}_train,{self.monitor_name}_val\n")
        self._data_files[estimator._computation_tree] = fname

    def on_fit_iter_end(
        self, *, estimator, node, from_reconstruction_attributes, data, **kwargs
    ):
        if node.depth != node.computation_tree.depth:
            return

        new_estimator = from_reconstruction_attributes

        X, y, X_val, y_val = data["X"], data["y"], data["X_val"], data["y_val"]

        if self.monitor == "objective_function":
            new_monitored_train, *_ = new_estimator.objective_function(
                X, y, normalize=True
            )
            if X_val is not None:
                new_monitored_val, *_ = new_estimator.objective_function(
                    X_val, y_val, normalize=True
                )
        elif callable(self.monitor):
            new_monitored_train = self.monitor(new_estimator, X, y)
            if X_val is not None:
                new_monitored_val = self.monitor(new_estimator, X_val, y_val)
        elif self.monitor is None or isinstance(self.monitor, str):
            from ..metrics import check_scoring

            scorer = check_scoring(new_estimator, self.monitor)
            new_monitored_train = scorer(new_estimator, X, y)
            if X_val is not None:
                new_monitored_val = scorer(new_estimator, X_val, y_val)

        if X_val is None:
            new_monitored_val = None

        with open(self._data_files[node.computation_tree], "a") as f:
            f.write(f"{node.idx},{new_monitored_train},{new_monitored_val}\n")

    def on_fit_end(self):
        pass

    # @property
    # def data(self):

    def plot(self):
        data_files = [p for p in Path(self._data_dir.name).iterdir() if p.is_file()]
        for f in data_files:
            data = pd.read_csv(f)
            fig, ax = plt.subplots()
            ax.plot(
                data["iteration"], data[f"{self.monitor_name}_train"], label="train set"
            )
            if self.on != "train_set":
                ax.plot(
                    data["iteration"],
                    data[f"{self.monitor_name}_val"],
                    label="validation set",
                )

            ax.set_xlabel("Number of iterations")
            ax.set_ylabel(self.monitor_name)

            ax.legend()
            plt.show()
