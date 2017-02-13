from __future__ import division, print_function

from copy import deepcopy

from .stats_node import StatsNode
from .criterion import impurity_improvement

FEATURE_THRESHOLD = 1e-7


class NewSplitter(object):
    """New type of splitter driven by data

    Parameters
    ----------
    X: ndarray, shape (n_samples, n_features)
        The dataset to be fitted.

    y: ndarray, shape (n_samples,)
        The residuals.

    sample_weight: ndarray, shape (n_samples,)
        The weight associated to each samples.

    sum_total_weighted_samples: float,
        The sum of the weights for all samples.

    split_record: SplitRecord,
        The split record to associate with this splitter.

    min_samples_leaf: int,
        The minimum number of samples to have a proper split record.

    min_weight_leaf: float,
        The minimum weight to have a proper split record.

    Attributes
    -------
    best_split_record: SplitRecord,
        The best split record found after iterating all the samples.
    """

    def __init__(self, X, y, sample_weight, sum_total_weighted_samples,
                 feature_idx, start_idx,
                 split_record,
                 min_samples_leaf, min_weight_leaf):
        # store the information related to the dataset
        self.X = X
        self.y = y
        self.sample_weight = sample_weight
        self.sum_total_weighted_samples = sum_total_weighted_samples

        # information about the feature and first sampled
        self.feature_idx = feature_idx
        self.start_idx = start_idx
        self.prev_idx = start_idx

        # split record to work with
        # make a deepcopy to not change the orignal object
        self.split_record = deepcopy(split_record)
        self.split_record.feature = self.feature_idx
        self.split_record.pos = self.start_idx
        # split to store the best split record
        self.best_split_record = deepcopy(split_record)

        # parameters for early stop of split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf

    def reset(self, feature_idx, start_idx, split_record):
        """Reset a splitter with a new samples set.
        It could correspond to a new feature to be scanned.
        """
        # information about the feature and first sampled
        self.feature_idx = feature_idx
        self.start_idx = start_idx
        self.prev_idx = start_idx

        # split record to work with
        # make a deepcopy to not change the original object
        self.split_record = deepcopy(split_record)
        self.split_record.feature = self.feature_idx
        self.split_record.pos = self.start_idx
        # split to store the best split record
        self.best_split_record = deepcopy(split_record)

    def update_stats(self, sample_idx):
        # make an update of the statistics
        # collect the statistics to add to the left node
        stats_samples = StatsNode(
            sum_residuals=self.y[sample_idx] * self.sample_weight[sample_idx],
            sum_sq_residuals=(self.y[sample_idx] ** 2.0 *
                              self.sample_weight[sample_idx]),
            n_samples=1,
            sum_weighted_samples=self.sample_weight[sample_idx])

        # add these statistics to the left child
        self.split_record.l_stats += stats_samples
        # update the statistics of the right child
        self.split_record.r_stats = (self.split_record.c_stats -
                                     self.split_record.l_stats)

    def node_evaluate_split(self, sample_idx):
        """Update the impurity and check the corresponding split should be
        kept.
        """
        feat_i = self.feature_idx

        self.update_stats(sample_idx)

        # check that the sample value are different enough
        change = abs(self.X[sample_idx, feat_i] -
                     self.X[self.split_record.pos, self.split_record.feature])
        if change > FEATURE_THRESHOLD:

            # check that there is enough samples to make the proper split
            if (self.split_record.l_stats.n_samples < self.min_samples_leaf or
                (self.split_record.r_stats.n_samples <
                 self.min_samples_leaf)):
                # early stopping
                return

            # check that the weight are enough important to make a proper split
            if ((self.split_record.l_stats.sum_weighted_samples <
                 self.min_weight_leaf) or
                (self.split_record.r_stats.sum_weighted_samples <
                 self.min_weight_leaf)):
                # early stopping
                return

            # compute the impurity improvement
            # FIXME we use the mse impurity for the moment
            c_impurity_improvement = impurity_improvement(
                self.split_record,
                self.sum_total_weighted_samples)

            # check the impurity improved
            if (c_impurity_improvement >
                    self.best_split_record.impurity_improvement):
                # reset the best split record
                threshold = (self.X[sample_idx, feat_i] +
                             self.X[self.prev_idx, feat_i]) / 2
                # update the best splitter
                self.best_split_record.reset(
                    feature=feat_i,
                    pos=sample_idx,
                    threshold=threshold,
                    impurity=self.split_record.impurity,
                    impurity_improvement=c_impurity_improvement,
                    nid=self.split_record.nid,
                    c_stats=self.split_record.c_stats,
                    l_stats=self.split_record.l_stats,
                    r_stats=self.split_record.r_stats)

        self.prev_idx = sample_idx
