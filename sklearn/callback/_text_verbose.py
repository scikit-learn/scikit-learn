# License: BSD 3 clause

import time

from . import BaseCallback


class TextVerbose(BaseCallback):

    auto_propagate = True
    request_stopping_criterion = True

    def __init__(self, min_time_between_calls=0):
        self.min_time_between_calls = min_time_between_calls

    def on_fit_begin(self, estimator, X=None, y=None):
        self.estimator = estimator
        self._start_time = time.perf_counter()

    def on_fit_iter_end(self, *, node, **kwargs):
        if node.depth != node.computation_tree.depth:
            return

        stopping_criterion = kwargs.get("stopping_criterion", None)
        tol = kwargs.get("tol", None)

        current_time = time.perf_counter() - self._start_time

        s = f"{node.description} {node.idx}"
        parent = node.parent
        while parent is not None and parent.parent is not None:
            s = f"{parent.description} {parent.idx} - {s}"
            parent = parent.parent

        msg = (
            f"[{parent.computation_tree.estimator_name}] {s} | time {current_time:.5f}s"
        )

        if stopping_criterion is not None and tol is not None:
            msg += f" | stopping_criterion={stopping_criterion:.3E} | tol={tol:.3E}"

        print(msg)

    def on_fit_end(self):
        pass
