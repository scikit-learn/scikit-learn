import numpy as np

_DEFAULT_TAGS = {
    'non_deterministic': False,
    'requires_positive_X': False,
    'requires_positive_y': False,
    'X_types': ['2darray'],
    'poor_score': False,
    'no_validation': False,
    'multioutput': False,
    "allow_nan": False,
    'stateless': False,
    'multilabel': False,
    '_skip_test': False,
    '_xfail_checks': False,
    'multioutput_only': False,
    'binary_only': False,
    'requires_fit': True,
    'preserves_dtype': [np.float64],
    'requires_y': False,
    'pairwise': False,
}


def _safe_tags(estimator, key=None):
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

    Returns
    -------
    tags : dict or tag value
        The estimator tags. A single value is returned if `key` is not None.
    """
    if hasattr(estimator, "_get_tags"):
        if key is not None:
            try:
                return estimator._get_tags().get(key, _DEFAULT_TAGS[key])
            except KeyError as exc:
                raise ValueError(
                    f"The key {key} is neither defined in _more_tags() in "
                    f"the class {repr(estimator)} nor a default tag key in "
                    f"_DEFAULT_TAGS."
                ) from exc
        else:
            return estimator._get_tags()
    else:
        if key is not None:
            try:
                default = _DEFAULT_TAGS[key]
            except KeyError as exc:
                raise ValueError(
                    f"The key {key} is not a default tags defined in "
                    f"_DEFAULT_TAGS and thus no default values are "
                    f"available. Use the parameter default if you want to "
                    f"define a default value."
                ) from exc
            return default
        else:
            return _DEFAULT_TAGS
