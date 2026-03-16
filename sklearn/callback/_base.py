# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from typing import Protocol, runtime_checkable


@runtime_checkable
class _BaseCallback(Protocol):
    """Protocol for the base callbacks."""

    def setup(self, context):
        """Method called at the beginning of the fit method of the estimator.

        For auto-propagated callbacks, this method is called only once, before running
        the fit method of the outermost estimator.

        Parameters
        ----------
        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task. This is usually the root context of the
            estimator but it can be an intermediate context if the estimator is a
            sub-estimator of a meta-estimator.
        """

    def teardown(self, context):
        """Method called after finishing the fit method of the estimator.

        For auto-propagated callbacks, this method is called only once, after finishing
        the fit method of the outermost estimator.

        Parameters
        ----------
        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task. This is usually the root context of the
            estimator but it can be an intermediate context if the estimator is a
            sub-estimator of a meta-estimator.
        """


@runtime_checkable
class FitCallback(_BaseCallback, Protocol):
    """Protocol for the callbacks evaluated on tasks during the fit of an estimator."""

    def on_fit_task_begin(self, context, **kwargs):
        """Method called at the beginning of each fit task of the estimator.

        Parameters
        ----------
        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task.

        **kwargs : dict
            Additional optional arguments holding information about the state of the
            fitting process at this task. The list of possible keys and corresponding
            values are described in detail at <TODO: add link>.
        """

    def on_fit_task_end(self, context, **kwargs):
        """Method called at the end of each fit task of the estimator.

        Parameters
        ----------
        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task.

        **kwargs : dict
            Additional optional arguments holding information about the state of the
            fitting process at this task. The list of possible keys and corresponding
            values are described in detail at <TODO: add link>.

        Returns
        -------
        stop : bool
            Whether or not to stop the current level of iterations at this task.
        """


@runtime_checkable
class AutoPropagatedCallback(_BaseCallback, Protocol):
    """Protocol for the auto-propagated callbacks

    An auto-propagated callback is a callback that is meant to be set on a top-level
    estimator and that is automatically propagated to its sub-estimators (if any).
    """

    @property
    def max_propagation_depth(self):
        """The maximum number of nested estimators at which the callback should be
        propagated.

        If set to None, the callback is propagated to sub-estimators at all nesting
        levels.
        """
