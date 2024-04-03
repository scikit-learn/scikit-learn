# License: BSD 3 clause
# Authors: the scikit-learn developers

from abc import ABC, abstractmethod


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
            Dictionary containing the training and validation data. The possible
            keys are "X_train", "y_train", "sample_weight_train", "X_val", "y_val"
            and "sample_weight_val".
        """

    @abstractmethod
    def on_fit_iter_end(self, estimator, task_node, **kwargs):
        """Method called at the end of each task of the estimator.

        Parameters
        ----------
        estimator : estimator instance
            The caller estimator. It might differ from the estimator passed to the
            `on_fit_begin` method for auto-propagated callbacks.

        task_node : TaskNode instance
            The caller task node.

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

    @abstractmethod
    def on_fit_end(self, task_node):
        """Method called at the end of the fit method of the estimator.

        Parameters
        ----------
        task_node : TaskNode instance
            The task node corresponding to the whole `fit` task. This is usually the
            root of the task tree of the estimator but it can be an intermediate node
            if the estimator is a sub-estimator of a meta-estimator.
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
