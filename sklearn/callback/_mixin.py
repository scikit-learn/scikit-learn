from ._base import CallbackProtocol
from ._callback_context import CallbackContext


class CallbackSupportMixin:
    """Mixin class to add callback support to an estimator."""

    def set_callbacks(self, callbacks):
        """Set callbacks for the estimator.

        Parameters
        ----------
        callbacks : callback or list of callbacks
            the callbacks to set.

        Returns
        -------
        self : estimator instance
            The estimator instance itself.
        """
        if not isinstance(callbacks, list):
            callbacks = [callbacks]

        if not all(isinstance(callback, CallbackProtocol) for callback in callbacks):
            raise TypeError("callbacks must follow the CallbackProtocol protocol.")

        self._skl_callbacks = callbacks

        return self

    def init_callback_context(self, task_name="fit"):
        """Initialize the callback context for the estimator.

        Parameters
        ----------
        task_name : str, default='fit'
            The name of the root task.

        Returns
        -------
        callback_fit_ctx : CallbackContext
            The callback context for the estimator.
        """
        # We don't initialize the callback context during _set_callbacks but in fit
        # because in the future we might want to have callbacks in predict/transform
        # which would require their own context.
        self._callback_fit_ctx = CallbackContext._from_estimator(
            estimator=self, task_name=task_name, task_id=0, max_tasks=1
        )

        return self._callback_fit_ctx
