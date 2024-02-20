# License: BSD 3 clause
# Authors: the scikit-learn developers

from abc import ABC, abstractmethod

# default values for the data dict passed to the callbacks
default_data = {
    "X_train": None,
    "y_train": None,
    "sample_weight_train": None,
    "X_val": None,
    "y_val": None,
    "sample_weight_val": None,
}


class BaseCallback(ABC):
    """Abstract class for the callbacks"""

    @abstractmethod
    def on_fit_begin(self, estimator, *, data):
        """Method called at the beginning of the fit method of the estimator.

        Parameters
        ----------
        estimator : estimator instance
            The estimator the callback is set on.

        data : dict
            Dictionary containing the training and validation data. The keys are
            "X_train", "y_train", "sample_weight_train", "X_val", "y_val",
            "sample_weight_val". The values are the corresponding data. If a key is
            missing, the corresponding value is None.
        """

    @abstractmethod
    def on_fit_end(self):
        """Method called at the end of the fit method of the estimator."""

    @abstractmethod
    def on_fit_iter_end(self, estimator, node, **kwargs):
        """Method called at the end of each computation node of the estimator.

        Parameters
        ----------
        estimator : estimator instance
            The caller estimator. It might differ from the estimator passed to the
            `on_fit_begin` method for auto-propagated callbacks.

        node : ComputationNode instance
            The caller computation node.

        **kwargs : dict
            arguments passed to the callback. Possible keys are

            - data: dict
                Dictionary containing the training and validation data. The keys are
                "X_train", "y_train", "sample_weight_train", "X_val", "y_val",
                "sample_weight_val". The values are the corresponding data. If a key is
                missing, the corresponding value is None.

            - stopping_criterion: float
                Usually iterations stop when `stopping_criterion <= tol`.
                This is only provided at the innermost level of iterations.

            - tol: float
                Tolerance for the stopping criterion.
                This is only provided at the innermost level of iterations.

            - from_reconstruction_attributes: estimator instance
                A ready to predict, transform, etc ... estimator as if the fit stopped
                at this node. Usually it's a copy of the caller estimator with the
                necessary attributes set but it can sometimes be an instance of another
                class (e.g. LogisticRegressionCV -> LogisticRegression)

            - fit_state: dict
                Model specific quantities updated during fit. This is not meant to be
                used by generic callbacks but by a callback designed for a specific
                estimator instead.

        Returns
        -------
        stop : bool or None
            Whether or not to stop the current level of iterations at this node.
        """

    @property
    def auto_propagate(self):
        """Whether or not this callback should be propagated to sub-estimators.

        An auto-propagated callback (from a meta-estimator to its sub-estimators) must
        be set on the meta-estimator. Its `on_fit_begin` and `on_fit_end` methods will
        only be called at the beginning and end of the fit method of the meta-estimator,
        while its `on_fit_iter_end` method will be called at each computation node of
        the meta-estimator and its sub-estimators.
        """
        return False

    def _is_propagated(self, estimator):
        """Check if this callback attached to estimator has been propagated from a
        meta-estimator.
        """
        return self.auto_propagate and hasattr(estimator, "_parent_node")

    # TODO: This is not used yet but will be necessary for next callbacks
    #       Uncomment when needed
    # @property
    # def request_stopping_criterion(self):
    #     return False

    # @property
    # def request_from_reconstruction_attributes(self):
    #     return False


class CallbackPropagatorMixin:
    """Mixin class for meta-estimators expected to propagate callbacks."""

    def _propagate_callbacks(self, sub_estimator, *, parent_node):
        """Propagate the auto-propagated callbacks to a sub-estimator.

        Parameters
        ----------
        sub_estimator : estimator instance
            The sub-estimator to propagate the callbacks to.

        parent_node : ComputationNode instance
            The computation node in this estimator to set as `parent_node` to the
            computation tree of the sub-estimator. It must be the node where the fit
            method of the sub-estimator is called.
        """
        if hasattr(sub_estimator, "_skl_callbacks") and any(
            callback.auto_propagate for callback in sub_estimator._skl_callbacks
        ):
            bad_callbacks = [
                callback.__class__.__name__
                for callback in sub_estimator._skl_callbacks
                if callback.auto_propagate
            ]
            raise TypeError(
                f"The sub-estimators ({sub_estimator.__class__.__name__}) of a"
                f" meta-estimator ({self.__class__.__name__}) can't have"
                f" auto-propagated callbacks ({bad_callbacks})."
                " Set them directly on the meta-estimator."
            )

        if not hasattr(self, "_skl_callbacks"):
            return

        propagated_callbacks = [
            callback for callback in self._skl_callbacks if callback.auto_propagate
        ]

        if not propagated_callbacks:
            return

        sub_estimator._parent_node = parent_node

        sub_estimator._set_callbacks(
            getattr(sub_estimator, "_skl_callbacks", []) + propagated_callbacks
        )


# Not a method of BaseEstimator because it might not be directly called from fit but
# by a non-method function called by fit
def _eval_callbacks_on_fit_iter_end(**kwargs):
    """Evaluate the `on_fit_iter_end` method of the callbacks.

    This function must be called at the end of each computation node.

    Parameters
    ----------
    kwargs : dict
        Arguments passed to the callback.

    Returns
    -------
    stop : bool
        Whether or not to stop the fit at this node.
    """
    estimator = kwargs.get("estimator")
    node = kwargs.get("node")

    if not hasattr(estimator, "_skl_callbacks") or node is None:
        return False

    # stopping_criterion and reconstruction_attributes can be costly to compute.
    # They are passed as lambdas for lazy evaluation. We only actually
    # compute them if a callback requests it.
    # TODO: This is not used yet but will be necessary for next callbacks
    #       Uncomment when needed
    # if any(cb.request_stopping_criterion for cb in estimator._skl_callbacks):
    #     kwarg = kwargs.pop("stopping_criterion", lambda: None)()
    #     kwargs["stopping_criterion"] = kwarg

    # if any(
    #         cb.request_from_reconstruction_attributes
    #         for cb in estimator._skl_callbacks
    # ):
    #     kwarg = kwargs.pop("from_reconstruction_attributes", lambda: None)()
    #     kwargs["from_reconstruction_attributes"] = kwarg

    return any(
        callback.on_fit_iter_end(**kwargs) for callback in estimator._skl_callbacks
    )
