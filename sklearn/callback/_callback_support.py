# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import copy

from sklearn.callback._base import Callback
from sklearn.callback._callback_context import CallbackContext


class CallbackSupportMixin:
    """Mixin class to add callback support to an estimator."""

    def set_callbacks(self, callbacks):
        """Set callbacks for the estimator.

        Parameters
        ----------
        callbacks : callback or list of callbacks
            The callbacks to set.

        Returns
        -------
        self : estimator instance
            The estimator instance itself.
        """
        if not isinstance(callbacks, list):
            callbacks = [callbacks]

        if not all(isinstance(callback, Callback) for callback in callbacks):
            raise TypeError("callbacks must follow the Callback protocol.")

        self._skl_callbacks = callbacks

        return self

    def _init_callback_context(self, task_name="fit", task_id=0, max_subtasks=0):
        """Initialize the callback context for the estimator.

        This method should be called once, at the beginning of the fit method.

        Parameters
        ----------
        task_name : str, default="fit"
            The name of the root task.

        task_id : int or str, default=0
            Identifier for the root task.

        max_subtasks : int or None, default=0
            The maximum number of subtasks that can be children of the root task. None
            means the maximum number of subtasks is not known in advance. 0 means it's a
            leaf.

        Returns
        -------
        callback_fit_ctx : CallbackContext
            The root callback context for the estimator.
        """
        self._callback_fit_ctx = CallbackContext._from_estimator(
            estimator=self,
            task_name=task_name,
            task_id=task_id,
            max_subtasks=max_subtasks,
        )

        return self._callback_fit_ctx

    def _from_reconstruction_attributes(self, *, reconstruction_attributes):
        """Return a copy of this estimator as if it was fitted.

        Parameters
        ----------
        reconstruction_attributes : callable
            A callable that has no arguments and returns the necessary fitted attributes
            to create a working fitted estimator from this instance.

            Using a callable allows lazy evaluation of the potentially costly
            reconstruction attributes.

        Returns
        -------
        fitted_estimator : estimator instance
            The fitted copy of this estimator.
        """
        new_estimator = copy.copy(self)  # XXX deepcopy ?
        for key, val in reconstruction_attributes().items():
            setattr(new_estimator, key, val)
        return new_estimator
