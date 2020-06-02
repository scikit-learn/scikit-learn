def _make_feature_names(n_features, prefix='x', input_features=None):
    """Make feature name strings from n_features.

    Either returns input_feature names if it is not None, or creates
    placeholder names based on n_features, by default,
    ['x0', 'x1', ..., 'xn_features'] is generated.

    Parameters
    ----------
    n_features : int
        Number of feature names to generate
    prefix : str, default='x'
        Prefix for each feature name.
    input_features : array-like of string
        Optional existing input features, returned unchanged if not None.

    Returns
    -------
    feature_names : list of str
        Generated feature names of length n_features.
    """
    if input_features is not None:
        return input_features
    return ["{}{}".format(prefix, i)
            for i in range(n_features)]
