# License: BSD 3 clause
# Authors: the scikit-learn developers

from sklearn.base import BaseEstimator, _fit_context, clone
from sklearn.callback import BaseCallback, CallbackPropagatorMixin
from sklearn.callback._base import _eval_callbacks_on_fit_iter_end
from sklearn.utils.parallel import Parallel, delayed


class TestingCallback(BaseCallback):
    def on_fit_begin(self, estimator, *, data):
        pass

    def on_fit_end(self):
        pass

    def on_fit_iter_end(self, estimator, node, **kwargs):
        pass


class TestingAutoPropagatedCallback(TestingCallback):
    auto_propagate = True


class NotValidCallback:
    """Unvalid callback since it does not inherit from `BaseCallback`."""

    def on_fit_begin(self, estimator, *, data):
        pass  # pragma: no cover

    def on_fit_end(self):
        pass  # pragma: no cover

    def on_fit_iter_end(self, estimator, node, **kwargs):
        pass  # pragma: no cover


class Estimator(BaseEstimator):
    _parameter_constraints: dict = {}

    def __init__(self, max_iter=20):
        self.max_iter = max_iter

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X, y):
        root = self._eval_callbacks_on_fit_begin(
            tree_structure=[
                {"stage": "fit", "n_children": self.max_iter},
                {"stage": "iter", "n_children": None},
            ],
            data={"X_train": X, "y_train": y},
        )

        for i in range(self.max_iter):
            if _eval_callbacks_on_fit_iter_end(
                estimator=self,
                node=root.children[i],
            ):
                break

        self.n_iter_ = i + 1

        return self


class MetaEstimator(BaseEstimator, CallbackPropagatorMixin):
    _parameter_constraints: dict = {}

    def __init__(
        self, estimator, n_outer=4, n_inner=3, n_jobs=None, prefer="processes"
    ):
        self.estimator = estimator
        self.n_outer = n_outer
        self.n_inner = n_inner
        self.n_jobs = n_jobs
        self.prefer = prefer

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X, y):
        root = self._eval_callbacks_on_fit_begin(
            tree_structure=[
                {"stage": "fit", "n_children": self.n_outer},
                {"stage": "outer", "n_children": self.n_inner},
                {"stage": "inner", "n_children": None},
            ],
            data={"X_train": X, "y_train": y},
        )

        Parallel(n_jobs=self.n_jobs, prefer=self.prefer)(
            delayed(_func)(self, self.estimator, X, y, node)
            for _, node in enumerate(root.children)
        )

        return self


def _func(meta_estimator, inner_estimator, X, y, parent_node):
    for _, node in enumerate(parent_node.children):
        est = clone(inner_estimator)
        meta_estimator._propagate_callbacks(est, parent_node=node)
        est.fit(X, y)

        _eval_callbacks_on_fit_iter_end(estimator=meta_estimator, node=node)

    _eval_callbacks_on_fit_iter_end(estimator=meta_estimator, node=parent_node)
