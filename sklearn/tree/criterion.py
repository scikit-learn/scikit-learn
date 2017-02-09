from __future__ import division, print_function


def _impurity_mse(stats_node):
    """Compute the impurity in MSE sense

    Parameters
    ----------
    stats_node: StatsNode,
        The split record from which to compute the impurity

    Returns
    -------
    impurity: float,
        The impurity.
    """
    impurity = (stats_node.sum_sq_residuals /
                stats_node.sum_weighted_samples)
    impurity -= ((stats_node.sum_residuals /
                  stats_node.sum_weighted_samples) ** 2.0)

    return impurity


def impurity_improvement(c_split_record, criterion='mse'):
    """Compute the impurity improvement.

    Parameters
    ----------
    p_split_record: SplitRecord,
        Split record of the parent node.

    c_split_record: SplitRecord,
        Split record of the current node.

    criterion: str, optional(default='mse')
        The criterion to use to compute the impurity improvement.

    Returns
    -------
    impurity_improvement: float,
        Impurity improvement
    """
    if criterion == 'mse':
        _impurity_func = _impurity_mse
    else:
        raise NotImplementedError

    # impurity current node
    c_impurity = _impurity_func(c_split_record.c_stats)
    # impurity left child
    l_impurity = _impurity_func(c_split_record.l_stats)
    # impurity right child
    r_impurity = _impurity_func(c_split_record.r_stats)

    # FIXME The formula is for the moment incorrect
    # we need to take into account the total number of samples to get a
    # relative improvement
    return c_impurity - l_impurity - r_impurity
