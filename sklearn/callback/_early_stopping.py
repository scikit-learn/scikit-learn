# License: BSD 3 clause

from urllib import request
from . import BaseCallback


class EarlyStopping(BaseCallback):

    request_from_reconstruction_attributes = True

    def __init__(
        self,
        X_val=None,
        y_val=None,
        monitor="objective_function",
        max_no_improvement=10,
        tol=1e-2,
    ):
        self.X_val = X_val
        self.y_val = y_val
        self.monitor = monitor
        self.max_no_improvement = max_no_improvement
        self.tol = tol

    def on_fit_begin(self, estimator, X=None, y=None):
        self.estimator = estimator
        self._no_improvement = {}
        self._last_monitored = {}

    def on_fit_iter_end(self, *, estimator, node, **kwargs):
        new_estimator = kwargs.get("from_reconstruction_attributes", None)

        if node.depth != self.estimator._computation_tree.depth:
            return

        if self.monitor == "objective_function":
            objective_function = kwargs.get("objective_function", None)
            monitored, *_ = objective_function(self.X_val)
        elif self.monitor == "TODO":
            pass

        if node.parent not in self._last_monitored or monitored < self._last_monitored[
            node.parent
        ] * (1 - self.tol):
            self._no_improvement[node.parent] = 0
            self._last_monitored[node.parent] = monitored
        else:
            self._no_improvement[node.parent] += 1

        if self._no_improvement[node.parent] >= self.max_no_improvement:
            return True

    def on_fit_end(self):
        pass
