# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import time

from sklearn.base import BaseEstimator, _fit_context, clone
from sklearn.callback import CallbackSupportMixin, with_fit_callbacks
from sklearn.callback._callback_support import get_callback_manager
from sklearn.utils.parallel import Parallel, delayed


class TestingCallback:
    """A minimal callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback doesn't define `max_propagation_depth` and is therefore not an
    `AutoPropagatedCallback`: it should not be propagated to sub-estimators.
    """

    def __init__(self):
        self.record = get_callback_manager().list()

    def setup(self, context):
        self.record.append({"name": "setup", "context": context})

    def on_fit_task_begin(self, context, **kwargs):
        self.record.append(
            {"name": "on_fit_task_begin", "context": context, "kwargs": kwargs}
        )

    def on_fit_task_end(self, context, **kwargs):
        self.record.append(
            {"name": "on_fit_task_end", "context": context, "kwargs": kwargs}
        )

    def teardown(self, context):
        self.record.append({"name": "teardown", "context": context})

    def count_hooks(self, hook_name):
        return len([rec for rec in self.record if rec["name"] == hook_name])


class TestingAutoPropagatedCallback(TestingCallback):
    """A minimal auto-propagated callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback defines `max_propagation_depth` and is therefore an
    `AutoPropagatedCallback`: it should be set on a top-level estimator and propagated
    to sub-estimators.
    """

    max_propagation_depth = None


class NotValidCallback:
    """Invalid callback since it's missing methods from the protocol."""

    def setup(self, estimator):
        pass  # pragma: no cover

    def on_fit_task_end(self, context, **kwargs):
        pass  # pragma: no cover


class FailingCallback(TestingCallback):
    """A callback that raises an error at some point."""

    def __init__(self, fail_at=None):
        super().__init__()
        self.fail_at = fail_at

    def setup(self, context):
        super().setup(context)
        if self.fail_at == "setup":
            raise ValueError("Failing callback failed at setup")

    def on_fit_task_begin(self, context, **kwargs):
        super().on_fit_task_begin(context, **kwargs)
        if self.fail_at == "on_fit_task_begin":
            raise ValueError("Failing callback failed at on_fit_task_begin")

    def on_fit_task_end(self, context, **kwargs):
        super().on_fit_task_end(context, **kwargs)
        if self.fail_at == "on_fit_task_end":
            raise ValueError("Failing callback failed at on_fit_task_end")

    def teardown(self, context):
        super().teardown(context)
        if self.fail_at == "teardown":
            raise ValueError("Failing callback failed at teardown")


class MaxIterEstimator(CallbackSupportMixin, BaseEstimator):
    """A class that mimics the behavior of an estimator.

    The iterative part uses a loop with a max number of iterations known in advance.
    """

    _parameter_constraints: dict = {}

    def __init__(self, max_iter=20, computation_intensity=0.001):
        self.max_iter = max_iter
        self.computation_intensity = computation_intensity

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X=None, y=None):
        callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
        callback_ctx.eval_on_fit_task_begin()

        for i in range(self.max_iter):
            subcontext = callback_ctx.subcontext(task_id=i).eval_on_fit_task_begin()

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.eval_on_fit_task_end(
                data={"X_train": X, "y_train": y},
            ):
                break

        callback_ctx.eval_on_fit_task_end(
            data={"X_train": X, "y_train": y},
        )

        self.n_iter_ = i + 1

        return self


class WhileEstimator(CallbackSupportMixin, BaseEstimator):
    """A class that mimics the behavior of an estimator.

    The iterative part uses a while loop with a number of iterations unknown in
    advance.
    """

    _parameter_constraints: dict = {}

    def __init__(self, computation_intensity=0.001):
        self.computation_intensity = computation_intensity

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X=None, y=None):
        callback_ctx = self._init_callback_context(max_subtasks=None)
        callback_ctx.eval_on_fit_task_begin()

        i = 0
        while True:
            subcontext = callback_ctx.subcontext(task_id=i).eval_on_fit_task_begin()

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.eval_on_fit_task_end(
                data={"X_train": X, "y_train": y},
            ):
                break

            if i == 20:
                break

            i += 1

        callback_ctx.eval_on_fit_task_end(
            data={"X_train": X, "y_train": y},
        )

        return self


class ThirdPartyEstimator(CallbackSupportMixin, BaseEstimator):
    """A class that mimics a third-party estimator with callback support only using
    public API.
    """

    def __init__(self, max_iter=20, computation_intensity=0.001):
        self.max_iter = max_iter
        self.computation_intensity = computation_intensity

    @with_fit_callbacks
    def fit(self, X=None, y=None):
        callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
        callback_ctx.eval_on_fit_task_begin()

        for i in range(self.max_iter):
            subcontext = callback_ctx.subcontext(task_id=i).eval_on_fit_task_begin()

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.eval_on_fit_task_end(
                data={"X_train": X, "y_train": y},
            ):
                break

        callback_ctx.eval_on_fit_task_end(
            data={"X_train": X, "y_train": y},
        )

        self.n_iter_ = i + 1

        return self


class ParentFitEstimator(MaxIterEstimator):
    """A class that mimics an estimator using its parent fit method."""

    _parameter_constraints: dict = {}

    def __init__(self, max_iter=20, computation_intensity=0.001):
        super().__init__(max_iter, computation_intensity)

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(self, X=None, y=None):
        return super().fit(X, y)


class NoCallbackEstimator(BaseEstimator):
    """A class that mimics an estimator without callback support."""

    def __init__(self, max_iter=20, computation_intensity=0.001):
        self.max_iter = max_iter
        self.computation_intensity = computation_intensity

    def fit(self, X=None, y=None):
        for i in range(self.max_iter):
            time.sleep(self.computation_intensity)  # Computation intensive task

        return self


class MetaEstimator(CallbackSupportMixin, BaseEstimator):
    """A class that mimics the behavior of a meta-estimator.

    It has two levels of iterations. The outer level uses parallelism and the inner
    level is done in a function that is not a method of the class. That function must
    therefore receive the estimator and the callback context as arguments.
    """

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
        callback_ctx = self._init_callback_context(max_subtasks=self.n_outer)
        callback_ctx.eval_on_fit_task_begin()

        Parallel(n_jobs=self.n_jobs, prefer=self.prefer)(
            delayed(_fit_subestimator)(
                self,
                self.estimator,
                X,
                y,
                outer_callback_ctx=callback_ctx.subcontext(
                    task_name="outer", task_id=i, max_subtasks=self.n_inner
                ),
            )
            for i in range(self.n_outer)
        )

        callback_ctx.eval_on_fit_task_end(
            data={"X_train": X, "y_train": y},
        )

        return self


def _fit_subestimator(meta_estimator, inner_estimator, X, y, *, outer_callback_ctx):
    outer_callback_ctx.eval_on_fit_task_begin()

    for i in range(meta_estimator.n_inner):
        est = clone(inner_estimator)

        inner_ctx = outer_callback_ctx.subcontext(task_name="inner", task_id=i)
        inner_ctx.propagate_callback_context(sub_estimator=est)
        inner_ctx.eval_on_fit_task_begin()

        est.fit(X, y)

        inner_ctx.eval_on_fit_task_end(
            data={"X_train": X, "y_train": y},
        )

    outer_callback_ctx.eval_on_fit_task_end(
        data={"X_train": X, "y_train": y},
    )


class NoSubtaskEstimator(CallbackSupportMixin, BaseEstimator):
    """A class mimicking an estimator without subtasks in fit."""

    @with_fit_callbacks
    def fit(self, X=None, y=None):
        callback_ctx = self._init_callback_context().eval_on_fit_task_begin()

        # No task performed

        callback_ctx.eval_on_fit_task_end()

        return self
