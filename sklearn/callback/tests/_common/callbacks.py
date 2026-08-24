# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Callbacks used to test the callback machinery."""

from sklearn.callback._transport import open_listener, send


class RecordingCallback:
    """A minimal callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback doesn't define `max_propagation_depth` and is therefore not an
    `AutoPropagatedCallback`: it should not be propagated to sub-estimators.
    """

    def __init__(self):
        self.record = []
        self._listener_handle = open_listener(self.record.append, owner=self)

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
        metadata=None,
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
                "estimator": estimator,
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
        send(
            self._listener_handle,
            {"name": "teardown", "estimator": estimator, "context": context},
        )

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
