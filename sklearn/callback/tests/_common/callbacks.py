# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Callbacks used to test the callback machinery."""

import copy

from sklearn.callback import AutoPropagatedCallback, FitCallback
from sklearn.callback._transport import open_listener, send


class RecordingCallback(FitCallback):
    """A minimal callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback doesn't define `max_propagation_depth` and is therefore not an
    `AutoPropagatedCallback`: it should not be propagated to sub-estimators.
    """

    def __init__(self):
        self.record = []
        self._listener_handle = open_listener(self.record.append, owner=self)

    def __getstate__(self):
        # We never access the record from a worker so there's no need to ship it.
        # Emptying it avoids nesting the record, callback and context in each other
        # which would grow the pickle size exponentially.
        return {**self.__dict__, "record": []}

    def setup(self, estimator, context):
        send(
            self._listener_handle,
            {"name": "setup", "estimator": estimator, "context": context},
        )

    def on_fit_task_begin(
        self,
        estimator,
        context,
        *,
        X=None,
        y=None,
        requested_arg_begin=None,
        fitted_estimator=None,
    ):
        send(
            self._listener_handle,
            {
                "name": "on_fit_task_begin",
                "estimator": estimator,
                "context": context,
                "kwargs": {
                    "X": X,
                    "y": y,
                    "requested_arg_begin": copy.copy(requested_arg_begin),
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
        requested_arg_end=None,
        fitted_estimator=None,
    ):
        send(
            self._listener_handle,
            {
                "name": "on_fit_task_end",
                "estimator": estimator,
                "context": context,
                "kwargs": {
                    "X": X,
                    "y": y,
                    "requested_arg_end": copy.copy(requested_arg_end),
                    "fitted_estimator": fitted_estimator,
                },
            },
        )

    def teardown(self, estimator, context):
        send(
            self._listener_handle,
            {"name": "teardown", "estimator": estimator, "context": context},
        )

    def count_hooks(self, hook_name):
        return len([rec for rec in self.record if rec["name"] == hook_name])


class RecordingAutoPropagatedCallback(RecordingCallback, AutoPropagatedCallback):
    """A minimal auto-propagated callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback defines `max_propagation_depth` and is therefore an
    `AutoPropagatedCallback`: it should be set on a top-level estimator and propagated
    to sub-estimators.
    """

    _max_propagation_depth = None


class NotValidCallback:
    """Invalid callback since it's not inheriting from FitCallback."""


class NotValidSetupPositionalCallback(RecordingCallback):
    """Invalid callback since it has invalid positional parameters."""

    def setup(self, estimator, context, not_valid_arg=None):
        pass  # pragma: no cover


class NotValidSetupKwargOnlyCallback(RecordingCallback):
    """Invalid callback since it has invalid kwarg-only parameters."""

    def setup(self, estimator, context, *, not_valid_kwarg=None):
        pass  # pragma: no cover


class NotValidFitTaskBeginCallback(RecordingCallback):
    """Invalid callback since it has invalid keyword-only parameters."""

    def on_fit_task_begin(self, estimator, context, not_valid_kwarg=None):
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


class SampleWeightCallback(RecordingCallback):
    """A callback that accepts sample_weight in its hooks."""

    def on_fit_task_begin(
        self, estimator, context, *, X=None, y=None, sample_weight=None
    ):
        super().on_fit_task_begin(
            estimator, context, X=X, y=y, requested_arg_begin=sample_weight
        )

    def on_fit_task_end(
        self, estimator, context, *, X=None, y=None, sample_weight=None
    ):
        super().on_fit_task_end(
            estimator, context, X=X, y=y, requested_arg_end=sample_weight
        )

    def _accept_sample_weight(self, hook_name):
        return hook_name in ("on_fit_task_begin", "on_fit_task_end")
