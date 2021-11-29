from ...base import is_classifier
from ...exceptions import NotFittedError
from ...utils.multiclass import type_of_target
from ...utils.validation import check_is_fitted


def _check_estimator_and_target_is_binary(estimator, y, target_type=None):
    """Helper to check that estimator is a binary classifier and y is binary.

    Parameters
    ----------
    estimator : estimator instance
        An estimator that should be used to predict the target.

    y : ndarray
        The associated target.

    target_type : str, default=None
        The type of the target `y` as returned by
        :func:`~sklearn.utils.multiclass.type_of_target`. If `None`, the type
        will be inferred by calling :func:`~sklearn.utils.multiclass.type_of_target`.
        Providing the type of the target could save time by avoid calling the
        :func:`~sklearn.utils.multiclass.type_of_target` function.

    Raises
    ------
    ValueError
        If the estimator or the target are not binary.
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

    if target_type is None:
        target_type = type_of_target(y)
    if target_type != "binary":
        raise ValueError(
            f"The target y is not binary. Got {target_type} type of target."
        )
