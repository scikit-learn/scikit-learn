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


# class NewSplitRecord(object):
#     """New split record
#     """

#     def __init__(self):
#         self.threshold = np.nan
#         self.feature = 0
#         self.pos = 0
#         self.impurity = np.inf

#         # statistic of the current node
#         self.sum_residual = 0
#         self.sq_sum_residual = 0
#         self.n_samples = 0
#         self.weighted_samples = 0

#         # statistic which will be used by the left and right child
#         # left node
#         self.sum_left_residual = 0
#         self.sq_sum_left_residual = 0
#         self.n_left_samples = 0
#         self.weighted_left_samples = 0
#         # right node
#         self.sum_right_residual = 0
#         self.sq_sum_right_residual = 0
#         self.n_right_samples = 0
#         self.weighted_right_samples = 0

# class NewCriterion(object):
#     """Criterion to compute the impurity improvement.
#     """

#     def __init__(self, y, sample_weight,
#                  p_split_record, c_split_record):

#         self.y = y
#         self.sample_weight = sample_weight

#         # contain the information of the parent split
#         self.p_split_record = p_split_record

#         # sum of residual stats
#         self.sum_left = 0
#         self.sum_right = 0
#         self.sum_total = self.p_split_record.sum_residual

#         self.sq_sum_left = 0
#         self.sq_sum_right = 0
#         self.sq_sum_total = self.p_split_record.sq_sum_residual

#         # sum the number of samples
#         self.n_total_samples = self.p_split_record.n_samples
#         self.n_left_samples = 0
#         self.n_right_samples = self.n_total_samples

#         # sum of weight stats
#         self.weighted_n_node_samples = self.p_split_record.weighted_samples
#         self.weighted_n_left = 0
#         self.weighted_n_right = 0

#         # needed during improvement impurity computation
#         self.weighted_n_total_samples = weighted_n_total_samples

#         # variable to keep the impurity
#         self.impurity_left = np.inf
#         self.impurity_right = np.inf

    # def update_stats(self, sample_idx):
    #     """ Update the statistics when a new sample is given.
    #     """
    #     # udpate the sum of the residual for the left node
    #     # sum_l += w_l * y_l
    #     self.sum_left += self.y[sample_idx] * self.sample_weight[sample_idx]
    #     # deduce the stat for the right node
    #     self.sum_right = self.sum_total - self.sum_left

    #     # update the sum of the residual squared for the left node
    #     # sq_sum_l += w_l * y_l * y_l
    #     self.sq_sum_left += (self.y[sample_idx] ** 2.0 *
    #                          self.sample_weight[sample_idx])
    #     # deduce the stat for the right node
    #     self.sq_sum_right = self.sq_sum_total - self.sq_sum_left

    #     # increment the number of samples in the left node
    #     self.n_left_samples += 1
    #     self.n_right_samples -= 1

    #     # compute the weights for the left node
    #     # sum_w_l += w_l
    #     self.weighted_n_left += self.sample_weight[sample_idx]
    #     # deduce the weights for the right node
    #     self.weighted_n_right = (self.weighted_n_node_samples -
    #                              self.weighted_n_left)

    # def impurity_improvement(self):
    #     """Compute the impurity improvement using the current statistics
    #     stored in the criterion.
    #     """
    #     # compute the left node impurity
    #     # i_l = sum of w_l * y_l * y_l - (sum w_l * y_l / ) ** 2
    #     self.impurity_left = (self.sq_sum_left / self.weighted_n_left -
    #                           (self.sum_left / self.weighted_n_left) ** 2.0)
    #     # compute the impurity for the right node
    #     self.impurity_right = (self.sq_sum_right / self.weighted_n_right -
    #                            (self.sum_right / self.weighted_n_right) ** 2.0)

    #     # compute the impurity improvement
    #     return (self.p_split_record.impurity - self.impurity_left -
    #             self.impurity_right)
