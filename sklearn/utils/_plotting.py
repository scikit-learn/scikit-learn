import numpy as np

from . import check_consistent_length, check_matplotlib_support
from .multiclass import type_of_target
from .validation import _check_pos_label_consistency
from ._response import _get_response_values_binary


class _BinaryClassifierCurveDisplayMixin:
    """Mixin class to be used in Displays requiring a binary classifier.

    The aim of this class is to centralize some validations regarding the estimator and
    the target and gather the response of the estimator.
    """

    def _validate_plot_params(self, *, ax=None, name=None):
        check_matplotlib_support(f"{self.__class__.__name__}.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()

        name = self.estimator_name if name is None else name
        return ax, ax.figure, name

    @classmethod
    def _validate_and_get_response_values(
        cls, estimator, X, y, *, response_method="auto", pos_label=None, name=None
    ):
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        name = estimator.__class__.__name__ if name is None else name

        y_pred, pos_label = _get_response_values_binary(
            estimator,
            X,
            response_method=response_method,
            pos_label=pos_label,
        )

        return y_pred, pos_label, name

    @classmethod
    def _validate_from_predictions_params(
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


def _validate_score_name(score_name, scoring, negate_score):
    """Validate the `score_name` parameter.

    If `score_name` is provided, we just return it as-is.
    Otherwise, we infer it from `scoring` by removing `"neg_"` if present and
    and make it English-readable. If `scoring` is `None`, then we return
    the default value `"Negative score"` if `negate_score` is `True` and
    `"Score"` otherwise.
    """
    if score_name is not None:
        return score_name
    elif scoring is None:
        return "Negative score" if negate_score else "Score"
    else:
        score_name = scoring.__name__ if callable(scoring) else scoring
        if negate_score and score_name.startswith("neg_"):
            score_name = score_name[4:]
        return score_name.replace("_", " ").capitalize()


def _compute_scale_type_ratio(data):
    """Compute the ratio between the largest and smallest inter-point distances.

    A value higher than 5 would indicate that the data was sampled from a
    log-uniform distribution.
    """
    diff = np.diff(np.sort(data))
    return diff.max() / diff.min()
