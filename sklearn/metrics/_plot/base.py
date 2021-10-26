from .._base import _check_pos_label_consistency
from .._classification import check_consistent_length

from ...base import is_classifier
from ...exceptions import NotFittedError
from ...utils import _get_response, check_matplotlib_support
from ...utils.multiclass import type_of_target
from ...utils.validation import check_is_fitted


class BaseBinaryClassifierCurveDisplay:
    """xxx"""

    def plot(self, *, name=None):
        check_matplotlib_support(f"{self.__class__.__name__}.plot")
        name = self.estimator_name if name is None else name
        return name

    @classmethod
    def from_estimator(
        cls, estimator, X, y, *, response_method="auto", pos_label=None, name=None
    ):
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        try:
            check_is_fitted(estimator)
        except NotFittedError as e:
            raise NotFittedError(
                f"This {estimator.__class__.__name__} instance is not fitted yet. Call "
                "'fit' with appropriate arguments before intending to use "
                f"'{cls.__name__}.from_estimator'."
            ) from e

        if not is_classifier(estimator):
            raise ValueError(
                f"This {cls.__name__} can only be used with a classifier. "
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

        if response_method == "auto":
            response_method = ["predict_proba", "decision_function"]

        name = estimator.__class__.__name__ if name is None else name

        y_pred, pos_label = _get_response(
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
