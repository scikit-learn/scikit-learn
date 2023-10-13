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
tag_sentinel = object()


def _safe_tags(estimator, key=None, default=tag_sentinel):
    """Safely get estimator tags.

    :class:`~sklearn.BaseEstimator` provides the estimator tags machinery.
    However, if an estimator does not inherit from this base class, we should
    fall-back to the default tags.

    For scikit-learn built-in estimators, we should still rely on
    `self._get_tags()`. `_safe_tags(est)` should be used when we are not sure
    where `est` comes from: typically `_safe_tags(self.base_estimator)` where
    `self` is a meta-estimator, or in the common checks.

    Parameters
    ----------
    estimator : estimator object
        The estimator from which to get the tag.

    key : str, default=None
        Tag name to get. By default (`None`), all tags are returned.

    default : obj, default=tag_sentinel
        Default value when the `key` tag is missing for `estimator`.

    Returns
    -------
    tags : dict or tag value
        The estimator tags. A single value is returned if `key` is not None.
    """
    default_is_sentinel = default is tag_sentinel

    if hasattr(estimator, "_get_tags"):
        tags_provider = "_get_tags()"
        tags = estimator._get_tags()
    elif hasattr(estimator, "_more_tags"):
        tags_provider = "_more_tags()"
        more_tags = estimator._more_tags()
        if default_is_sentinel:
            tags = {**_DEFAULT_TAGS, **more_tags}
        else:
            tags = more_tags
    else:
        tags_provider = "_DEFAULT_TAGS"
        tags = _DEFAULT_TAGS if default_is_sentinel else {}

    if key is not None:
        if key not in tags:
            if default_is_sentinel:
                raise ValueError(
                    f"The key {key} is not defined in {tags_provider} for the "
                    f"class {estimator.__class__.__name__}."
                )
            return default
        return tags[key]
    return tags
