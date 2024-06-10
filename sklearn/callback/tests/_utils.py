# License: BSD 3 clause
# Authors: the scikit-learn developers

import time

from sklearn.base import BaseEstimator, _fit_context, clone
from sklearn.callback import CallbackSupportMixin
from sklearn.utils.parallel import Parallel, delayed


class TestingCallback:
    def _on_fit_begin(self, estimator, *, data):
        pass

    def _on_fit_end(self):
        pass

    def _on_fit_iter_end(self, estimator, node, **kwargs):
        pass


class TestingAutoPropagatedCallback(TestingCallback):
    max_estimator_depth = None


class NotValidCallback:
    """Unvalid callback since it's missing a method from the protocol.'"""

    def _on_fit_begin(self, estimator, *, data):
        pass  # pragma: no cover

    def _on_fit_iter_end(self, estimator, node, **kwargs):
        pass  # pragma: no cover


class Estimator(CallbackSupportMixin, BaseEstimator):
    _parameter_constraints: dict = {}

    def __init__(self, max_iter=20, computation_intensity=0.001):
        self.max_iter = max_iter
        self.computation_intensity = computation_intensity

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X=None, y=None):
        callback_ctx = self.init_callback_context().eval_on_fit_begin(
            estimator=self, data={"X_train": X, "y_train": y}
        )

        for i in range(self.max_iter):
            subcontext = callback_ctx.subcontext(task_id=i, max_tasks=self.max_iter)

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.eval_on_fit_iter_end(
                estimator=self,
                data={"X_train": X, "y_train": y},
            ):
                break

        self.n_iter_ = i + 1

        return self


class WhileEstimator(CallbackSupportMixin, BaseEstimator):
    _parameter_constraints: dict = {}

    def __init__(self, computation_intensity=0.001):
        self.computation_intensity = computation_intensity

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X=None, y=None):
        callback_ctx = self.init_callback_context().eval_on_fit_begin(
            estimator=self, data={"X_train": X, "y_train": y}
        )

        i = 0
        while True:
            subcontext = callback_ctx.subcontext(task_id=i, max_tasks=None)

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.eval_on_fit_iter_end(
                estimator=self,
                data={"X_train": X, "y_train": y},
            ):
                break

            if i == 20:
                break

            i += 1

        return self


class MetaEstimator(CallbackSupportMixin, BaseEstimator):
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
    def fit(self, X=None, y=None):
        callback_ctx = self.init_callback_context().eval_on_fit_begin(
            estimator=self, data={"X_train": X, "y_train": y}
        )

        Parallel(n_jobs=self.n_jobs, prefer=self.prefer)(
            delayed(_func)(
                self,
                self.estimator,
                X,
                y,
                callback_ctx=callback_ctx.subcontext(
                    task_name="outer", task_id=i, max_tasks=self.n_outer
                ),
            )
            for i in range(self.n_outer)
        )

        return self


def _func(meta_estimator, inner_estimator, X, y, *, callback_ctx):
    for i in range(meta_estimator.n_inner):
        est = clone(inner_estimator)

        inner_ctx = callback_ctx.subcontext(
            task_name="inner", task_id=i, max_tasks=meta_estimator.n_inner
        ).propagate_callbacks(sub_estimator=est)

        est.fit(X, y)

        inner_ctx.eval_on_fit_iter_end(
            estimator=meta_estimator,
            data={"X_train": X, "y_train": y},
        )

    callback_ctx.eval_on_fit_iter_end(
        estimator=meta_estimator,
        data={"X_train": X, "y_train": y},
    )
