# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from typing import Protocol, runtime_checkable


@runtime_checkable
class Callback(Protocol):
    """Protocol for the callbacks."""

    def on_fit_begin(self, estimator):
        """Method called at the beginning of the fit method of the estimator.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.
        """

    def on_fit_task_end(self, estimator, context, **kwargs):
        """Method called at the end of each task of the estimator.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook. It might differ from the estimator
            passed to the `_on_fit_begin` method for auto-propagated callbacks.

        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task.

        **kwargs : dict
            Information about the state of the fitting process at this task. Possible
            keys are:

            - "data": dict
                Dictionary containing the training and validation data. The possible
                keys are "X_train", "y_train", "sample_weight_train", "X_val", "y_val",
                "sample_weight_val".

            - "stopping_criterion": float
                Usually iterations stop when `stopping_criterion <= tol`.
                This is only provided at the innermost level of iterations.

            - "tol": float
                Tolerance for the stopping criterion.
                This is only provided at the innermost level of iterations.

            - "from_reconstruction_attributes": estimator instance
                A ready to predict, transform, etc ... estimator as if the fit stopped
                at the end of this task.

            - "fit_state": dict
                Model specific quantities updated during fit. This is not meant to be
                used by generic callbacks but by a callback designed for a specific
                estimator instead.

        Returns
        -------
        stop : bool
            Whether or not to stop the current level of iterations at this task.
        """

    def on_fit_end(self, estimator, context):
        """Method called at the end of the fit method of the estimator.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.

        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task. This is usually the root context of the
            estimator but it can be an intermediate context if the estimator is a
            sub-estimator of a meta-estimator.
        """


@runtime_checkable
class AutoPropagatedCallback(Protocol):
    """Protocol for the auto-propagated callbacks

    An auto-propagated callback is a callback that is meant to be set on a top-level
    estimator and that is automatically propagated to its sub-estimators (if any).
    """

    @property
    def max_estimator_depth(self):
        """The maximum number of nested estimators at which the callback should be
        propagated.

        If set to None, the callback is propagated to sub-estimators at all nesting
        levels.
        """
