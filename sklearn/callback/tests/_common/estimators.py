# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Estimators used to test the callback machinery."""

import time

import numpy as np

from sklearn.base import BaseEstimator, _fit_context, clone
from sklearn.callback import CallbackSupportMixin, with_callbacks
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _manual_routing,
    _routing_enabled,
    process_routing,
)
from sklearn.utils.parallel import Parallel, delayed


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
    def fit(self, X=None, y=None, sample_weight=None, **fit_params):
        callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
        if _routing_enabled():
            routed_params = process_routing(
                self, "fit", sample_weight=sample_weight, **fit_params
            )
        else:
            routed_params = _manual_routing(
                {"callback": self._get_manual_callback_params(sample_weight)}
            )
        callback_ctx.call_on_fit_task_begin(
            estimator=self,
            X=X,
            y=y,
            metadata=routed_params.callback.on_fit_task_begin,
        )

        for i in range(self.max_iter):
            subcontext = callback_ctx.subcontext(task_name=f"iteration {i}")
            subcontext.call_on_fit_task_begin(
                estimator=self,
                X=X,
                y=y,
                metadata=routed_params.callback.on_fit_task_begin,
            )

            time.sleep(self.computation_intensity)  # Computation intensive task

            if subcontext.call_on_fit_task_end(
                estimator=self,
                X=X,
                y=y,
                reconstruction_attributes=lambda: {"n_iter_": i + 1},
                metadata=routed_params.callback.on_fit_task_end,
            ):
                break

        self.n_iter_ = i + 1

        callback_ctx.call_on_fit_task_end(
            estimator=self,
            X=X,
            y=y,
            reconstruction_attributes={},
            metadata=routed_params.callback.on_fit_task_end,
        )

        return self

    def predict(self, X):
        return np.mean(X, axis=1) * self.n_iter_

    def get_metadata_routing(self):
        router = MetadataRouter(owner=self).add_self_request(self)
        return self._add_callback_routing(router)


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

            if subcontext.call_on_fit_task_end(
                estimator=self,
                X=X,
                y=y,
                reconstruction_attributes={"n_iter_": i + 1},
            ):
                break

            if i == 20:
                break

            i += 1

        self.n_iter_ = i + 1

        callback_ctx.call_on_fit_task_end(
            estimator=self, X=X, y=y, reconstruction_attributes={}
        )

        return self

    def predict(self, X):
        return np.mean(X, axis=1) * self.n_iter_


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

    def predict(self, X):
        return np.zeros(X.shape[0])


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
    def fit(self, X=None, y=None, **fit_params):
        callback_ctx = self._init_callback_context(
            max_subtasks=self.n_outer, sequential_subtasks=False
        )
        routed_params = process_routing(self, "fit", **fit_params)
        callback_ctx.call_on_fit_task_begin(
            estimator=self,
            X=X,
            y=y,
            metadata=routed_params.callback.on_fit_task_begin,
        )

        outer_callback_contexts = [
            callback_ctx.subcontext(
                task_name="outer", task_id=i, max_subtasks=self.n_inner
            )
            for i in range(self.n_outer)
        ]

        Parallel(n_jobs=self.n_jobs, prefer=self.prefer)(
            delayed(_fit_subestimator)(
                self,
                self.estimator,
                X=X,
                y=y,
                fit_params=routed_params.estimator.fit,
                callback_params=routed_params.callback,
                outer_callback_ctx=outer_callback_contexts[i],
            )
            for i in range(self.n_outer)
        )

        callback_ctx.call_on_fit_task_end(
            estimator=self,
            X=X,
            y=y,
            metadata=routed_params.callback.on_fit_task_end,
        )

        return self

    def get_metadata_routing(self):
        router = MetadataRouter(owner=self).add(
            estimator=self.estimator,
            method_mapping=MethodMapping().add(caller="fit", callee="fit"),
        )
        return self._add_callback_routing(router)


def _fit_subestimator(
    meta_estimator,
    inner_estimator,
    *,
    X,
    y,
    fit_params,
    callback_params,
    outer_callback_ctx,
):
    outer_callback_ctx.call_on_fit_task_begin(
        estimator=meta_estimator,
        X=X,
        y=y,
        metadata=callback_params.on_fit_task_begin,
    )

    for i in range(meta_estimator.n_inner):
        est = clone(inner_estimator)

        inner_ctx = outer_callback_ctx.subcontext(task_name="inner")
        with inner_ctx.propagate_callback_context(est):
            inner_ctx.call_on_fit_task_begin(
                estimator=meta_estimator,
                X=X,
                y=y,
                metadata=callback_params.on_fit_task_begin,
            )

            est.fit(X=X, y=y, **fit_params)

            inner_ctx.call_on_fit_task_end(
                estimator=meta_estimator,
                X=X,
                y=y,
                metadata=callback_params.on_fit_task_end,
            )

    outer_callback_ctx.call_on_fit_task_end(
        estimator=meta_estimator,
        X=X,
        y=y,
        metadata=callback_params.on_fit_task_end,
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
