from .._base import _check_pos_label_consistency
from .._classification import check_consistent_length

from ...base import is_classifier
from ...exceptions import NotFittedError
from ...utils import _get_response_values, check_matplotlib_support
from ...utils.multiclass import type_of_target
from ...utils.validation import check_is_fitted


def _check_estimator_target(estimator, y):
    """Helper to check that estimator is a binary classifier and y is binary.

    This function is aside from the class `BinaryClassifierCurveDisplayMixin`
    below because it allows to have consistent error messages between the
    displays and the plotting functions.

    FIXME: Move into `BinaryClassifierCurveDisplayMixin.from_estimator` when
    the plotting functions will be removed in 1.2.
    """
    try:
        check_is_fitted(estimator)
    except NotFittedError as e:
        raise NotFittedError(
            f"This {estimator.__class__.__name__} instance is not fitted yet. Call "
            "'fit' with appropriate arguments before intending to use it to plotting "
            "functionalities."
        ) from e

    if not is_classifier(estimator):
        raise ValueError(
            "This plotting functionalities only support a binary classifier. "
            f"Got a {estimator.__class__.__name__} instead."
        )
    elif len(estimator.classes_) != 2:
        raise ValueError(
            f"This {estimator.__class__.__name__} instance is not a binary "
            "classifier. It was fitted on multiclass problem with "
            f"{len(estimator.classes_)} classes."
        )
    elif type_of_target(y) != "binary":
        raise ValueError(
            f"The target y is not binary. Got {type_of_target(y)} type of target."
        )


class BinaryClassifierCurveDisplayMixin:
    """Mixin class to be used in Displays requiring a binary classifier.

    The aim of this class is to make some validations regarding the estimator and the
    target and gather the response of the estimator.
    """

    def plot(self, *, name=None):
        """Validate the curve label.

        This method will check that matplotlib is installed and validate
        the name to be shown in the legend of the plot.

        Parameters
        ----------
        name : str, default=None
            Name introduced in the legend of the plot to point out the current
            estimator.

        Returns
        -------
        name : str
            The validated name to be shown in the legend.

        Raises
        ------
        ImportError
            If matplotlib is not installed.
        """
        check_matplotlib_support(f"{self.__class__.__name__}.plot")
        name = self.estimator_name if name is None else name
        return name

    @classmethod
    def from_estimator(
        cls, estimator, X, y, *, response_method="auto", pos_label=None, name=None
    ):
        """Validate estimator and target and gather response.

        This method will check that matplotlib is installed. Then, it will
        validate that `estimator` is a binary classifier and `y` is binary.
        Then, the response values will be gather from the estimator after
        checking that the response method is available in the estimator.

        Parameters
        ----------
        estimator : estimator instance
            Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a classifier.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Target values.

        response_method : {'predict_proba', 'decision_function', 'auto'} \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response. If set to 'auto',
            :term:`predict_proba` is tried first and if it does not exist
            :term:`decision_function` is tried next.

        pos_label : str or int, default=None
            The class considered as the positive class when computing the roc
            auc metrics. By default, `estimators.classes_[1]` is considered as
            the positive class.

        name : str, default=None
            Name to label the curve. If `None`, use the name of the
            estimator.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,)
            The response values given by the estimator.

        pos_label : str, int or None
            The class considered as the positive class when computing
            the metrics. Returns `None` if `estimator` is a regressor.

        name : str
            Name to label the curve.

        Raises
        ------
        ImportError
            If matplotlib is not installed.
        """
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        _check_estimator_target(estimator, y)

        if response_method == "auto":
            response_method = ["predict_proba", "decision_function"]

        name = estimator.__class__.__name__ if name is None else name

        y_pred, pos_label = _get_response_values(
            estimator,
            X,
            y,
            response_method,
            pos_label=pos_label,
        )

        return y_pred, pos_label, name

    @classmethod
    def from_predictions(
        cls, y_true, y_pred, *, sample_weight=None, pos_label=None, name=None
    ):
        """Validate true and predicted targets and sample weight.

        The length of the true and predicted targets as well as sample weights
        will be validated.

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_pred : array-like of shape (n_samples,)
            Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
            (as returned by “decision_function” on some classifiers).

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        pos_label : str or int, default=None
            The label of the positive class. When `pos_label=None`, if `y_true`
            is in {-1, 1} or {0, 1}, `pos_label` is set to 1, otherwise an
            error will be raised.

        name : str, default=None
            Name to label the curve. If `None`, name will be set to
            `"Classifier"`.

        Returns
        -------
        pos_label : int
            If `pos_label` can be inferred, it will be returned.

        name : str
            Name to label the curve.

        Raises
        ------
        ImportError
            If matplotlib is not installed.
        """
        check_matplotlib_support(f"{cls.__name__}.from_predictions")

        if type_of_target(y_true) != "binary":
            raise ValueError(
                f"The target y is not binary. Got {type_of_target(y_true)} type of"
                " target."
            )

        check_consistent_length(y_true, y_pred, sample_weight)
        pos_label = _check_pos_label_consistency(pos_label, y_true)

        name = name if name is not None else "Classifier"

        return pos_label, name
