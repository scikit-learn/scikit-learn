import inspect
import warnings

import numpy as np

_DEFAULT_TAGS = {
    "array_api_support": False,
    "non_deterministic": False,
    "requires_positive_X": False,
    "requires_positive_y": False,
    "X_types": ["2darray"],
    "poor_score": False,
    "no_validation": False,
    "multioutput": False,
    "allow_nan": False,
    "stateless": False,
    "multilabel": False,
    "_skip_test": False,
    "_xfail_checks": False,
    "multioutput_only": False,
    "binary_only": False,
    "requires_fit": True,
    "preserves_dtype": [np.float64],
    "requires_y": False,
    "pairwise": False,
}


# TODO(1.7) Remove `_more_tags` support
def _walk_mro_more_tags(estimator):
    collected_tags = {}
    for base_class in reversed(inspect.getmro(estimator.__class__)):
        if hasattr(base_class, "_more_tags"):
            # need the if because mixins might not have _more_tags
            # but might do redundant work in estimators
            # (i.e. calling more tags on BaseEstimator multiple times)
            more_tags = base_class._more_tags(estimator)
            collected_tags.update(more_tags)
    return collected_tags


def _safe_tags(estimator, key=None):
    """Safely get estimator tags.

    :class:`~sklearn.BaseEstimator` provides the estimator tags machinery.
    However, if an estimator does not inherit from this base class, we should
    fall-back to the default tags.

    For scikit-learn built-in estimators, we should still rely on
    `self.__sklearn_tags__()`. `_safe_tags(est)` should be used when we are not sure
    where `est` comes from: typically `_safe_tags(self.base_estimator)` where
    `self` is a meta-estimator, or in the common checks.

    Parameters
    ----------
    estimator : estimator object
        The estimator from which to get the tag.

    key : str, default=None
        Tag name to get. By default (`None`), all tags are returned.

    Returns
    -------
    tags : dict or tag value
        The estimator tags. A single value is returned if `key` is not None.
    """
    if hasattr(estimator, "__sklearn_tags__"):
        tags_provider = "__sklearn_tags__()"
        tags = estimator.__sklearn_tags__()
    elif hasattr(estimator, "_get_tags"):
        # TODO(1.7) Remove `_get_tags` support
        warnings.warn(
            "_get_tags() was deprecated in 1.1 support will be removed in 1.3. "
            "Please use __sklearn_tags__ instead.",
            FutureWarning,
        )
        tags_provider = "_get_tags()"
        tags = estimator._get_tags()
    elif hasattr(estimator, "_more_tags"):
        # TODO(1.7) Remove `_more_tags` support
        warnings.warn(
            "_more_tags() was deprecated in 1.1 support will be removed in 1.3. "
            "Please use __sklearn_tags__ instead.",
            FutureWarning,
        )
        tags_provider = "_more_tags()"
        tags = {**_DEFAULT_TAGS, **_walk_mro_more_tags(estimator)}
    else:
        tags_provider = "_DEFAULT_TAGS"
        tags = _DEFAULT_TAGS

    if key is not None:
        if key not in tags:
            raise ValueError(
                f"The key {key} is not defined in {tags_provider} for the "
                f"class {estimator.__class__.__name__}."
            )
        return tags[key]
    return tags
