# License: BSD 3 clause

from . import BaseCallback


class EarlyStopping(BaseCallback):
    request_from_reconstruction_attributes = True

    def __init__(
        self,
        monitor="objective_function",
        on="validation_set",
        higher_is_better=False,
        validation_split="auto",
        max_no_improvement=10,
        threshold=1e-2,
    ):
        from ..model_selection import KFold

        self.validation_split = validation_split
        if validation_split == "auto":
            self.validation_split = KFold(n_splits=5, shuffle=True, random_state=42)
        self.monitor = monitor
        self.on = on
        self.higher_is_better = higher_is_better
        self.max_no_improvement = max_no_improvement
        self.threshold = threshold

        self._no_improvement = {}
        self._last_monitored = {}

    def on_fit_begin(self, estimator, X=None, y=None):
        pass

    def on_fit_iter_end(self, *, estimator, node, **kwargs):
        if node.depth != node.computation_tree.depth:
            return

        reconstructed_estimator = kwargs.pop("from_reconstruction_attributes")
        data = kwargs.pop("data")

        X = data["X_val"] if self.on == "validation_set" else data["X"]
        y = data["y_val"] if self.on == "validation_set" else data["y"]

        if self.monitor == "objective_function":
            new_monitored, *_ = reconstructed_estimator.objective_function(
                X, y, normalize=True
            )
        elif callable(self.monitor):
            new_monitored = self.monitor(reconstructed_estimator, X, y)
        elif self.monitor is None or isinstance(self.monitor, str):
            from ..metrics import check_scoring

            scorer = check_scoring(reconstructed_estimator, self.monitor)
            new_monitored = scorer(reconstructed_estimator, X, y)

        if self._score_improved(node, new_monitored):
            self._no_improvement[node.parent] = 0
            self._last_monitored[node.parent] = new_monitored
        else:
            self._no_improvement[node.parent] += 1

        if self._no_improvement[node.parent] >= self.max_no_improvement:
            return True

    def _score_improved(self, node, new_monitored):
        if node.parent not in self._last_monitored:
            return True

        last_monitored = self._last_monitored[node.parent]
        if self.higher_is_better:
            return new_monitored > last_monitored * (1 + self.threshold)
        else:
            return new_monitored < last_monitored * (1 - self.threshold)

    def on_fit_end(self):
        pass

    @property
    def request_validation_split(self):
        return self.on == "validation_set"
