def _get_feature_index(fx, feature_names=None):
    """Get feature index.

    Parameters
    ----------
    fx : int, str or bool
        Feature index, name or mask.

    pos : int
        Position of the feature in the feature list.

    feature_names : list of str, default=None
        All feature names from which to search the indices.

    Returns
    -------
    idx : int
        Feature index.
    """
    if isinstance(fx, str):
        if feature_names is None:
            raise ValueError(
                "When the feature is a string, `feature_names` should be a "
                "list of feature names."
            )
        try:
            fx = feature_names.index(fx)
        except ValueError as e:
            raise ValueError(f"Feature {fx} not in feature_names") from e
    return int(fx)
