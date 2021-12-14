# License: BSD 3 clause

from copy import copy
from pathlib import Path
from tempfile import mkdtemp
import time

import matplotlib.pyplot as plt
import pandas as pd

from . import BaseCallback


class ConvergenceMonitor(BaseCallback):
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

    request_reconstruction_attributes = True

    def __init__(self, *, monitor="objective_function", X_val=None, y_val=None):
        self.X_val = X_val
        self.y_val = y_val
        self._data_file = Path(mkdtemp()) / "convergence_monitor.csv"

    def on_fit_begin(self, estimator, *, X=None, y=None):
        self.estimator = estimator
        self.X_train = X
        self.y_train = y
        self._start_time = {}

    def on_fit_iter_end(self, *, node, **kwargs):
        if node.depth != node.computation_tree.depth:
            return

        reconstruction_attributes = kwargs.get("reconstruction_attributes", None)
        if reconstruction_attributes is None:
            return

        new_estimator = copy(self.estimator)
        for key, val in reconstruction_attributes.items():
            setattr(new_estimator, key, val)

        if node.idx == 0:
            self._start_time[node.parent] = time.perf_counter()
            curr_time = 0
        else:
            curr_time = time.perf_counter() - self._start_time[node.parent]

        obj_train, *_ = new_estimator.objective_function(self.X_train, self.y_train, normalize=True)
        if self.X_val is not None:
            obj_val, *_ = new_estimator.objective_function(self.X_val, self.y_val, normalize=True)
        else:
            obj_val = None

        ancestors = node.get_ancestors()[:0:-1]
        ancestors_desc = [
            f"{n.computation_tree.estimator_name}-{n.description}" for n in ancestors
        ]
        ancestors_idx = [f"{n.idx}" for n in ancestors]

        if not self._data_file.exists():
            with open(self._data_file, "w") as f:
                f.write(
                    f"{','.join(ancestors_desc)},iteration,time,obj_train,obj_val\n"
                )

        with open(self._data_file, "a") as f:
            f.write(
                f"{','.join(ancestors_idx)},{node.idx},{curr_time},{obj_train},{obj_val}\n"
            )

    def on_fit_end(self):
        pass

    def get_data(self):
        if not hasattr(self, "data"):
            self.data = pd.read_csv(self._data_file)
        return self.data

    def plot(self, x="iteration"):
        data = self.get_data()

        # all columns but iteration, time, obj_train, obj_val
        group_by_columns = list(data.columns[:-4])
        groups = data.groupby(group_by_columns)

        for key in groups.groups.keys():
            group = groups.get_group(key)
            fig, ax = plt.subplots()

            ax.plot(group[x], group["obj_train"], label="obj_train")
            if self.X_val is not None:
                ax.plot(group[x], group["obj_val"], label="obj_val")

            if x == "iteration":
                x_label = "Number of iterations"
            elif x == "time":
                x_label = "Time (s)"
            ax.set_xlabel(x_label)
            ax.set_ylabel("objective function")

            ax.legend()
            plt.show()
