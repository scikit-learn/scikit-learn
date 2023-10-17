from functools import partial

from joblib.parallel import Parallel, delayed

from sklearn.base import BaseEstimator, _fit_context, clone
from sklearn.callback import BaseCallback
from sklearn.callback._base import _eval_callbacks_on_fit_iter_end


class TestingCallback(BaseCallback):
    def on_fit_begin(self, estimator, *, X=None, y=None):
        pass

    def on_fit_end(self):
        pass

    def on_fit_iter_end(self, estimator, node, **kwargs):
        pass


class TestingAutoPropagatedCallback(TestingCallback):
    auto_propagate = True


class NotValidCallback:
    def on_fit_begin(self, estimator, *, X=None, y=None):
        pass

    def on_fit_end(self):
        pass

    def on_fit_iter_end(self, estimator, node, **kwargs):
        pass


class Estimator(BaseEstimator):
    _parameter_constraints = {}

    def __init__(self, max_iter=20):
        self.max_iter = max_iter

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X, y):
        root, X, y, X_val, y_val = self._eval_callbacks_on_fit_begin(
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
                from_reconstruction_attributes=partial(
                    self._from_reconstruction_attributes,
                    reconstruction_attributes=lambda: {"n_iter_": i + 1},
                ),
                data={"X": X, "y": y, "X_val": X_val, "y_val": y_val},
            ):
                break

        self.n_iter_ = i + 1

        return self

    def objective_function(self, X, y=None, normalize=False):
        return 0, 0, 0


class MetaEstimator(BaseEstimator):
    _parameter_constraints = {}

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
        root, X, y, _, _ = self._eval_callbacks_on_fit_begin(
            levels=[
                {"descr": "fit", "max_iter": self.n_outer},
                {"descr": "outer", "max_iter": self.n_inner},
                {"descr": "inner", "max_iter": None},
            ],
            X=X,
            y=y,
        )

        Parallel(n_jobs=self.n_jobs, prefer=self.prefer)(
            delayed(self._func)(self.estimator, X, y, node, i)
            for i, node in enumerate(root.children)
        )

        return self

    def _func(self, estimator, X, y, parent_node, i):
        for j, node in enumerate(parent_node.children):
            est = clone(estimator)
            self._propagate_callbacks(est, parent_node=node)
            est.fit(X, y)

            _eval_callbacks_on_fit_iter_end(estimator=self, node=node)

        _eval_callbacks_on_fit_iter_end(estimator=self, node=parent_node)

        return
