# License: BSD 3 clause

from abc import ABC, abstractmethod


# Not a method of BaseEstimator because it might be called from an extern function
def _eval_callbacks_on_fit_iter_end(**kwargs):
    """Evaluate the on_fit_iter_end method of the callbacks

    This function should be called at the end of each computation node.

    Parameters
    ----------
    kwargs : dict
        arguments passed to the callback.

    Returns
    -------
    stop : bool
        Whether or not to stop the fit at this node.
    """
    estimator = kwargs.get("estimator")
    node = kwargs.get("node")

    if not hasattr(estimator, "_callbacks") or node is None:
        return False

    estimator._computation_tree._tree_status[node.tree_status_idx] = True

    # stopping_criterion and reconstruction_attributes can be costly to compute. They
    # are passed as lambdas for lazy evaluation. We only actually compute them if a
    # callback requests it.
    if any(
        getattr(callback, "request_stopping_criterion", False)
        for callback in estimator._callbacks
    ):
        kwarg = kwargs.pop("stopping_criterion", lambda: None)()
        kwargs["stopping_criterion"] = kwarg

    if any(
        getattr(callback, "request_reconstruction_attributes", False)
        for callback in estimator._callbacks
    ):
        kwarg = kwargs.pop("reconstruction_attributes", lambda: None)()
        kwargs["reconstruction_attributes"] = kwarg

    return any(callback.on_fit_iter_end(**kwargs) for callback in estimator._callbacks)


class BaseCallback(ABC):
    """Abstract class for the callbacks"""

    @abstractmethod
    def on_fit_begin(self, estimator, *, X=None, y=None):
        """Method called at the beginning of the fit method of the estimator

        Parameters
        ----------
        estimator: estimator instance
            The estimator the callback is set on.
        X: ndarray or sparse matrix, default=None
            The training data.
        y: ndarray, default=None
            The target.
        """
        pass

    @abstractmethod
    def on_fit_end(self):
        """Method called at the end of the fit method of the estimator"""
        pass

    @abstractmethod
    def on_fit_iter_end(self, estimator, node, **kwargs):
        """Method called at the end of each computation node of the estimator

        Parameters
        ----------
        estimator : estimator instance
            The caller estimator. It might differ from the estimator passed to the
            `on_fit_begin` method for auto-propagated callbacks.

        node : ComputationNode instance
            The caller computation node.

        kwargs : dict
            arguments passed to the callback. Possible keys are

            - stopping_criterion: float
                Usually iterations stop when `stopping_criterion <= tol`.
                This is only provided at the innermost level of iterations.

            - tol: float
                Tolerance for the stopping criterion.
                This is only provided at the innermost level of iterations.

            - reconstruction_attributes: dict
                Necessary attributes to construct an estimator (by copying this
                estimator and setting these as attributes) which will behave as if
                the fit stopped at this node.
                This is only provided at the outermost level of iterations.

            - fit_state: dict
                Model specific quantities updated during fit. This is not meant to be
                used by generic callbacks but by a callback designed for a specific
                estimator instead.

        Returns
        -------
        stop : bool or None
            Whether or not to stop the current level of iterations at this node.
        """
        pass


class AutoPropagatedMixin:
    """Mixin for auto-propagated callbacks

    An auto-propagated callback (from a meta-estimator to its sub-estimators) must be
    set on the meta-estimator. Its `on_fit_begin` and `on_fit_end` methods will only be
    called at the beginning and end of the fit method of the meta-estimator, while its
    `on_fit_iter_end` method will be called at each computation node of the
    meta-estimator and its sub-estimators.
    """

    pass
