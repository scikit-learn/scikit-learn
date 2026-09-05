# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._metadata_requests import _MetadataRequester


class BaseCallback:
    """Base class for the callbacks."""

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

    def _accept_sample_weight(self, hook_name):
        """Method to check if the instance accepts sample_weight for a given hook.

        Callbacks that can use sample_weight must overwrite this method to signal if
        a hook can be forwarded sample_weight when metadata routing is disabled.
        """
        # TODO(slep006): remove when metadata routing is always enabled.
        return False


class FitCallback(BaseCallback, _MetadataRequester):
    """Base class for the callbacks evaluated during the fit of an estimator."""

    def on_fit_task_begin(
        self,
        estimator,
        context,
        *,
        X=None,
        y=None,
        fitted_estimator=None,
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

        fitted_estimator : estimator instance
            A new instance of the estimator that is ready to predict, transform, etc ...
            as if fit had stopped at the beginning of this task.
        """

    def on_fit_task_end(
        self,
        estimator,
        context,
        *,
        X=None,
        y=None,
        fitted_estimator=None,
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

        fitted_estimator : estimator instance
            A new instance of the estimator that is ready to predict, transform, etc ...
            as if fit had stopped at the end of this task.

        Returns
        -------
        stop : bool
            Whether or not to stop the current level of iterations at this task.
        """


class AutoPropagatedCallback(BaseCallback):
    """Base class for the auto-propagated callbacks

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
        return self._max_propagation_depth

    @max_propagation_depth.setter
    def max_propagation_depth(self, value):
        """The setter: This is where you can add logic/validation"""
        if not isinstance(value, int) or value < 0:
            raise ValueError("max_propagation_depth must be a positive integer.")
        self._max_propagation_depth = value
