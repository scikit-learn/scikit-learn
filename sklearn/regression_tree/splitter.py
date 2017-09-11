from __future__ import division, print_function

from .stats_node import StatsNode
from .split_record import SplitRecord
from .criterion import impurity_improvement

FEATURE_THRESHOLD = 1e-7


class Splitter(object):
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
        self.split_record = SplitRecord()

        split_record.copy_to(self.split_record)
        self.split_record.feature = self.feature_idx
        self.split_record.pos = self.start_idx

        self.best_split_record = SplitRecord()
        # split to store the best split record
        split_record.copy_to(self.best_split_record)
        # set the right stats node equal to the left stats node
        split_record.c_stats.copy_to(split_record.r_stats)

        # parameters for early stop of split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf

        # create a stats_node to get information about a single sample
        # used in update_stats
        self.stats_samples = StatsNode(0., 0., 0, 0.)

    def reset(self, feature_idx, start_idx, split_record):
        """Reset a splitter with a new samples set.
        It could correspond to a new feature to be scanned.
        """
        # information about the feature and first sampled
        self.feature_idx = feature_idx
        self.start_idx = start_idx
        self.prev_idx = start_idx

        # split record to work with
        split_record.copy_to(self.split_record)
        self.split_record.feature = self.feature_idx
        self.split_record.pos = self.start_idx
        # split to store the best split record
        split_record.copy_to(self.best_split_record)
        # set the right stats node equal to the left stats node
        split_record.c_stats.copy_to(split_record.r_stats)

    def update_stats(self, sample_idx):

        # make an update of the statistics
        # collect the statistics to add to the left node
        self.stats_samples.reset(
            sum_y=self.y[sample_idx] * self.sample_weight[sample_idx],
            sum_sq_y=(self.y[sample_idx] ** 2.0 *
                      self.sample_weight[sample_idx]),
            n_samples=1,
            sum_weighted_samples=self.sample_weight[sample_idx])

        # add these statistics to the left child
        self.split_record.l_stats += self.stats_samples
        # update the statistics of the right child
        self.split_record.r_stats -= self.stats_samples

    def node_evaluate_split(self, sample_idx):
        """Update the impurity and check the corresponding split should be
        kept.
        """
        feat_i = self.feature_idx

        # check that the two consecutive samples are not the same
        b_samples_var = (abs(self.X[sample_idx, feat_i] -
                             self.X[self.prev_idx, feat_i]) >
                         FEATURE_THRESHOLD)

        # check that there is enough samples to make a split
        b_n_samples = not (
            self.split_record.l_stats.n_samples <
            self.min_samples_leaf or
            self.split_record.r_stats.n_samples <
            self.min_samples_leaf)

        # check that the weights corresponding to samples is great enough
        b_weight_samples = not(
            self.split_record.l_stats.sum_weighted_samples <
            self.min_weight_leaf or
            self.split_record.r_stats.sum_weighted_samples <
            self.min_weight_leaf)

        # try to split if necessary
        if b_samples_var and b_n_samples and b_weight_samples:

            # compute the impurity improvement
            # FIXME we use the mse impurity for the moment
            c_impurity_improvement = impurity_improvement(
                self.split_record.c_stats,
                self.split_record.l_stats,
                self.split_record.r_stats,
                self.sum_total_weighted_samples,
                'mse')

            # check the impurity improved
            if (c_impurity_improvement >
                    self.best_split_record.impurity_improvement):
                # update the best split
                self.best_split_record.reset(
                    feature=feat_i,
                    pos=self.prev_idx,
                    threshold=((self.X[sample_idx, feat_i] +
                                self.X[self.prev_idx, feat_i]) / 2.),
                    impurity=self.split_record.impurity,
                    impurity_improvement=c_impurity_improvement,
                    nid=self.split_record.nid,
                    c_stats=self.split_record.c_stats,
                    l_stats=self.split_record.l_stats,
                    r_stats=self.split_record.r_stats)

        self.update_stats(sample_idx)
        self.prev_idx = sample_idx
