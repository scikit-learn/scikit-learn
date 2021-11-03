from ...base import is_classifier
from ...exceptions import NotFittedError
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
