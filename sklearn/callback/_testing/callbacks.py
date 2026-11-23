# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Callbacks used to test the callback machinery."""

import copy

from sklearn.callback import AutoPropagatedCallbackMixin, FitCallback
from sklearn.callback._transport import Channel


class RecordingCallback(FitCallback):
    """A minimal callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback doesn't define `max_propagation_depth` and is therefore not an
    auto-propagated callback: it should not be propagated to sub-estimators.
    """

    def __init__(self):
        self.record = []
        self._channel = Channel(self.record.append)

    def __getstate__(self):
        # We never access the record from a worker so there's no need to ship it.
        # Emptying it avoids nesting the record, callback and context in each other
        # which would grow the pickle size exponentially.
        return {**self.__dict__, "record": []}

    def setup(self, estimator, context):
        self._channel.send(
            {"name": "setup", "estimator": estimator, "context": context},
        )

    def on_fit_task_begin(
        self,
        estimator,
        context,
        *,
        X=None,
        y=None,
        fitted_estimator=None,
        requested_arg_begin=None,
    ):
        """Method called at the beginning of each fit task of the estimator, which
        records the arguments passed to the hook.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.

        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task.

        X : array-like
            The training data at this task.

        y : array-like
            The training target values at this task.

        fitted_estimator : estimator instance
            A new instance of the estimator that is ready to predict, transform, etc ...
            as if fit had stopped at the beginning of this task.

        requested_arg_begin : object
            An argument that can be requested by the callback for this hook.
        """
        self._channel.send(
            {
                "name": "on_fit_task_begin",
                "estimator": estimator,
                "context": context,
                "kwargs": {
                    "X": X,
                    "y": y,
                    "fitted_estimator": fitted_estimator,
                    "requested_arg_begin": copy.copy(requested_arg_begin),
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
        fitted_estimator=None,
        requested_arg_end=None,
    ):
        """Method called at the beginning of each fit task of the estimator, which
        records the arguments passed to the hook.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.

        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task.

        X : array-like
            The training data at this task.

        y : array-like
            The training target values at this task.

        fitted_estimator : estimator instance
            A new instance of the estimator that is ready to predict, transform, etc ...
            as if fit had stopped at the beginning of this task.

        requested_arg_end : object
            An argument that can be requested by the callback for this hook.
        """
        self._channel.send(
            {
                "name": "on_fit_task_end",
                "estimator": estimator,
                "context": context,
                "kwargs": {
                    "X": X,
                    "y": y,
                    "fitted_estimator": fitted_estimator,
                    "requested_arg_end": copy.copy(requested_arg_end),
                },
            },
        )

    def teardown(self, estimator, context):
        self._channel.send(
            {"name": "teardown", "estimator": estimator, "context": context},
        )
        self._channel.disconnect()

    def count_hooks(self, hook_name):
        """Count the recorded calls for a given hook.

        Parameters
        ----------
        hook_name : str
            Name of the hook to count.

        Returns
        -------
        count : int
            Number of recorded calls for `hook_name`.
        """
        return len([rec for rec in self.record if rec["name"] == hook_name])


class RecordingAutoPropagatedCallback(AutoPropagatedCallbackMixin, RecordingCallback):
    """A minimal auto-propagated callback used for smoke testing purposes.

    This callback keeps a record of the hooks called for introspection.

    This callback defines `max_propagation_depth` and is therefore an
    auto-propagated callback: it should be set on a top-level estimator and propagated
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
