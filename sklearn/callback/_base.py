# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from typing import Protocol, runtime_checkable


@runtime_checkable
class _BaseCallback(Protocol):
    """Protocol for the base callbacks."""

    def setup(self, estimator, context):
        """Method called at the beginning of the fit method of the estimator.

        For auto-propagated callbacks, this method is called only once, before running
        the fit method of the outermost estimator.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.

        context : `sklearn.callback.CallbackContext` instance
            Context of the corresponding task. This is usually the root context of the
            estimator but it can be an intermediate context if the estimator is a
            sub-estimator of a meta-estimator.
        """

    def teardown(self, estimator, context):
        """Method called after finishing the fit method of the estimator.

        For auto-propagated callbacks, this method is called only once, after finishing
        the fit method of the outermost estimator.

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
class FitCallback(_BaseCallback, Protocol):
    """Protocol for the callbacks evaluated on tasks during the fit of an estimator."""

    def on_fit_task_begin(
        self, estimator, context, X=None, y=None, metadata=None, fitted_estimator=None
    ):
        """Method called at the beginning of each fit task of the estimator.

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

        metadata : dict
            Training metadata at this task, e.g. sample weights.

        fitted_estimator : estimator instance
            A new instance of the estimator that is ready to predict, transform, etc ...
            as if fit had stopped at the beginning of this task.
        """

    def on_fit_task_end(
        self, estimator, context, X=None, y=None, metadata=None, fitted_estimator=None
    ):
        """Method called at the end of each fit task of the estimator.

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

        metadata : dict
            Training metadata at this task, e.g. sample weights.

        fitted_estimator : estimator instance
            A new instance of the estimator that is ready to predict, transform, etc ...
            as if fit had stopped at the end of this task.

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
