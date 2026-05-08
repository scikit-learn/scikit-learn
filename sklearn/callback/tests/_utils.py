# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import time

import numpy as np

from sklearn.base import BaseEstimator, _fit_context, clone
from sklearn.callback import CallbackSupportMixin, with_callbacks
from sklearn.callback._transport import open_listener, send
from sklearn.utils.parallel import Parallel, delayed


class RecordingCallback:
    """A minimal callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback doesn't define `max_propagation_depth` and is therefore not an
    `AutoPropagatedCallback`: it should not be propagated to sub-estimators.
    """

    def __init__(self):
        self.record = []
        self._listener_handle = open_listener(self.record.append)

    def setup(self, estimator, context):
        send(self._listener_handle, {"name": "setup", "context": context})

    def on_fit_task_begin(
        self,
        estimator,
        context,
        *,
        X=None,
        y=None,
        metadata=None,
        fitted_estimator=None,
    ):
        send(
            self._listener_handle,
            {
                "name": "on_fit_task_begin",
                "context": context,
                "kwargs": {
                    "X": X,
                    "y": y,
                    "metadata": metadata,
                    "fitted_estimator": fitted_estimator,
                },
            },
        )

    def on_fit_task_end(
        self,
        estimator,
        context,
        *,
        X=None,
        y=None,
        metadata=None,
        fitted_estimator=None,
    ):
        send(
            self._listener_handle,
            {
                "name": "on_fit_task_end",
                "context": context,
                "kwargs": {
                    "X": X,
                    "y": y,
                    "metadata": metadata,
                    "fitted_estimator": fitted_estimator,
                },
            },
        )

    def teardown(self, estimator, context):
        send(self._listener_handle, {"name": "teardown", "context": context})

    def count_hooks(self, hook_name):
        return len([rec for rec in self.record if rec["name"] == hook_name])


class RecordingAutoPropagatedCallback(RecordingCallback):
    """A minimal auto-propagated callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback defines `max_propagation_depth` and is therefore an
    `AutoPropagatedCallback`: it should be set on a top-level estimator and propagated
    to sub-estimators.
    """

    max_propagation_depth = None


class NotValidCallback:
    """Invalid callback since it's missing methods from the protocol."""

    def setup(self, estimator, context):
        pass  # pragma: no cover

    def on_fit_task_end(self, estimator, context):
        pass  # pragma: no cover


class NotValidHookCallback(RecordingCallback):
    """Invalid callback since it has invalid parameters in the hooks signatures."""

    def on_fit_task_begin(self, estimator, context, *, not_valid_kwarg=None):
        pass  # pragma: no cover


class FailingCallback(RecordingCallback):
    """A callback that raises an error at some point."""

    def __init__(self, fail_at=None):
        super().__init__()
        self.fail_at = fail_at

    def setup(self, estimator, context):
        super().setup(estimator, context)
        if self.fail_at == "setup":
            raise ValueError("Failing callback failed at setup")

    def on_fit_task_begin(self, estimator, context):
        super().on_fit_task_begin(estimator, context)
        if self.fail_at == "on_fit_task_begin":
            raise ValueError("Failing callback failed at on_fit_task_begin")

    def on_fit_task_end(self, estimator, context):
        super().on_fit_task_end(estimator, context)
        if self.fail_at == "on_fit_task_end":
            raise ValueError("Failing callback failed at on_fit_task_end")

    def teardown(self, estimator, context):
        super().teardown(estimator, context)
        if self.fail_at == "teardown":
            raise ValueError("Failing callback failed at teardown")


class StopFitCallback(RecordingCallback):
    """A callback with a `on_fit_task_end` hook returning True."""

    def on_fit_task_end(self, estimator, context):
        super().on_fit_task_end(estimator, context)
        return True


class NotRequiredKwargsCallback(RecordingCallback):
    """A callback with a `on_fit_task_end` not requiring all possible kwargs."""

    def on_fit_task_end(self, estimator, context, *, X=None, y=None):
        super().on_fit_task_end(estimator, context, X=X, y=y)


class MaxIterEstimator(CallbackSupportMixin, BaseEstimator):
    """A class that mimics the behavior of an estimator.

    The iterative part uses a loop with a max number of iterations known in advance.

    This estimator computes arbitrary predictions by averaging the feature
    values and multiplying the result by the number of iterations done
    in fit.
    """

    _parameter_constraints: dict = {}

    def __init__(self, max_iter=20, computation_intensity=0.001):
        self.max_iter = max_iter
        self.computation_intensity = computation_intensity

    @_fit_context(prefer_skip_nested_validation=False)
    def fit(
        self,
        X=None,
        y=None,
        *,
        sample_weight=None,
    ):
        callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
        metadata = {"sample_weight": sample_weight} if sample_weight is not None else {}
        callback_ctx.call_on_fit_task_begin(estimator=self, X=X, y=y, metadata=metadata)

        for i in range(self.max_iter):
            subcontext = callback_ctx.subcontext(task_name=f"iteration {i}")
            subcontext.call_on_fit_task_begin(
                estimator=self, X=X, y=y, metadata=metadata
            )

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.call_on_fit_task_end(
                estimator=self,
                X=X,
                y=y,
                metadata=metadata,
                reconstruction_attributes=lambda: {"n_iter_": i + 1},
            ):
                break

        self.n_iter_ = i + 1

        callback_ctx.call_on_fit_task_end(
            estimator=self,
            X=X,
            y=y,
            metadata=metadata,
            reconstruction_attributes={},
        )

        return self

    def predict(self, X):
        return np.mean(X, axis=1) * self.n_iter_


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
        callback_ctx.call_on_fit_task_begin(estimator=self, X=X, y=y)

        i = 0
        while True:
            subcontext = callback_ctx.subcontext()
            subcontext.call_on_fit_task_begin(estimator=self, X=X, y=y)

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.call_on_fit_task_end(estimator=self, X=X, y=y):
                break

            if i == 20:
                break

            i += 1

        callback_ctx.call_on_fit_task_end(estimator=self, X=X, y=y)

        return self


class ThirdPartyEstimator(CallbackSupportMixin):
    """A class that mimics a third-party estimator with callback support only using
    public API.
    """

    def __init__(self, max_iter=20, computation_intensity=0.001):
        self.max_iter = max_iter
        self.computation_intensity = computation_intensity

    @with_callbacks
    def fit(self, X=None, y=None):
        callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
        callback_ctx.call_on_fit_task_begin(estimator=self, X=X, y=y)

        for i in range(self.max_iter):
            subcontext = callback_ctx.subcontext()
            subcontext.call_on_fit_task_begin(estimator=self, X=X, y=y)

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.call_on_fit_task_end(estimator=self, X=X, y=y):
                break

        callback_ctx.call_on_fit_task_end(estimator=self, X=X, y=y)

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
    def fit(self, X=None, y=None, sample_weight=None):
        callback_ctx = self._init_callback_context(
            max_subtasks=self.n_outer, sequential_subtasks=False
        )
        metadata = {"sample_weight": sample_weight} if sample_weight is not None else {}
        callback_ctx.call_on_fit_task_begin(estimator=self, X=X, y=y, metadata=metadata)

        Parallel(n_jobs=self.n_jobs, prefer=self.prefer)(
            delayed(_fit_subestimator)(
                self,
                self.estimator,
                X=X,
                y=y,
                metadata=metadata,
                outer_callback_ctx=callback_ctx.subcontext(
                    task_name="outer", task_id=i, max_subtasks=self.n_inner
                ),
            )
            for i in range(self.n_outer)
        )

        callback_ctx.call_on_fit_task_end(estimator=self, X=X, y=y, metadata=metadata)

        return self


def _fit_subestimator(
    meta_estimator, inner_estimator, *, X, y, metadata, outer_callback_ctx
):
    outer_callback_ctx.call_on_fit_task_begin(
        estimator=meta_estimator, X=X, y=y, metadata=metadata
    )

    for i in range(meta_estimator.n_inner):
        est = clone(inner_estimator)

        inner_ctx = outer_callback_ctx.subcontext(task_name="inner")
        with inner_ctx.propagate_callback_context(est):
            inner_ctx.call_on_fit_task_begin(
                estimator=meta_estimator, X=X, y=y, metadata=metadata
            )

            est.fit(X=X, y=y, **metadata)

            inner_ctx.call_on_fit_task_end(
                estimator=meta_estimator, X=X, y=y, metadata=metadata
            )

    outer_callback_ctx.call_on_fit_task_end(
        estimator=meta_estimator, X=X, y=y, metadata=metadata
    )


class HeterogeneousMetaEstimator(CallbackSupportMixin):
    """A meta-estimator that fits a list of estimators in order."""

    def __init__(self, estimators):
        self.estimators = estimators

    @with_callbacks
    def fit(self, X=None, y=None):
        callback_ctx = self._init_callback_context(max_subtasks=len(self.estimators))
        callback_ctx.call_on_fit_task_begin(estimator=self, X=X, y=y)

        for i, est in enumerate(self.estimators):
            task_name = f"fit {est.__class__.__name__}" if est else f"skip {i}"
            subcontext = callback_ctx.subcontext(task_name=task_name)
            if est is not None:
                est = clone(est)
                with subcontext.propagate_callback_context(est):
                    subcontext.call_on_fit_task_begin(estimator=self, X=X, y=y)
                    est.fit(X, y)
                    subcontext.call_on_fit_task_end(estimator=self, X=X, y=y)
            else:
                subcontext.call_on_fit_task_begin(estimator=self, X=X, y=y)
                subcontext.call_on_fit_task_end(estimator=self, X=X, y=y)

        callback_ctx.call_on_fit_task_end(estimator=self, X=X, y=y)

        return self


class NoSubtaskEstimator(CallbackSupportMixin, BaseEstimator):
    """A class mimicking an estimator without subtasks in fit."""

    @with_callbacks
    def fit(self, X=None, y=None):
        callback_ctx = self._init_callback_context().call_on_fit_task_begin(
            estimator=self, X=X, y=y
        )

        # No task performed

        callback_ctx.call_on_fit_task_end(estimator=self, X=X, y=y)

        return self
