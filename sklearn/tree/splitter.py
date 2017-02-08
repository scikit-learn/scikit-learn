from __future__ import division, print_function
from copy import copy

from .criterion import NewCriterion
from .criterion import NewSplitRecord

FEATURE_THRESHOLD = 1e-7

class NewSplitter(object):
    """New type of splitter driven by data
    """

    def __init__(self, X, y, sample_weight, feature_idx, start,
                 weighted_n_total_samples,
                 p_split_record, min_samples_leaf, min_weight_leaf):
        self.X = X
        self.y = y
        self.sample_weight = sample_weight
        self.feature_idx = feature_idx
        self.weighted_n_total_samples = weighted_n_total_samples

        # get the parent criterion
        self.p_split_record = copy(p_split_record)
        print("A new splitter n_samples is ", self.p_split_record.n_samples)

        # parameters for early stop of split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf

        # split record for the current node
        self.split_record = NewSplitRecord()
        self.split_record.feature = feature_idx
        self.split_record.pos = start
        self.split_record.n_samples = p_split_record.n_samples

        # create the criterion for the current splitter
        self.criterion = NewCriterion(self.X, self.y, self.sample_weight,
                                      self.weighted_n_total_samples,
                                      self.p_split_record)

    def reset(self, feature_idx, start, p_split_record):
        """Reset a splitter with a new samples set.
        It could correspond to a new feature to be scanned.
        """
        self.feature_idx = feature_idx
        self.p_split_record = copy(p_split_record)

        # reinitialize the split record
        self.split_record.__init__()
        self.split_record.feature = feature_idx
        self.split_record.pos = start
        #self.split_record.n_samples = p_split_record.n_samples
        #self.split_record.weighted_samples = p_split_record.n_samples

        # reinitialize the criterion
        self.criterion.__init__(self.X, self.y, self.sample_weight,
                                self.weighted_n_total_samples, p_split_record)

    def node_evaluate_split(self, sample_idx):
        """Update the impurity and check the corresponding split should be
        kept.
        """
        feat_i = self.feature_idx

        print('cur_sample/pos', self.X[sample_idx, feat_i],
              self.X[self.split_record.pos, self.split_record.feature])

        # make an update of the statistics
        self.criterion.update_stats(sample_idx)

        # check that the sample value are different enough
        if self.X[sample_idx, feat_i] != self.X[self.split_record.pos,
                                                self.split_record.feature]:


            # check that there is enough samples to make the proper split
            if (self.criterion.n_left_samples < self.min_samples_leaf or
                    self.criterion.n_right_samples < self.min_samples_leaf):
                # early stopping
                print("Stopping here. n_left or n_right less than min_samples_leaf")
                print("n_left/n_right - ", self.criterion.n_left_samples,
                      self.criterion.n_right_samples)
                return

            # check that the weight are enough important to make a proper split
            if (self.criterion.weighted_n_left < self.min_weight_leaf or
                    self.criterion.weighted_n_right < self.min_weight_leaf):
                # early stopping
                print("Stopping here. n_weighted_left or n_weighted_right "
                      "less than min_samples_leaf")
                print("n_w_left/n_w_right - ", self.criterion.weighted_n_left,
                      self.criterion.weighted_n_right)
                return

            # compute the impurity improvement
            current_impurity = self.criterion.impurity_improvement()
            print('current_imp spl_rec_imp',
                  current_impurity, self.split_record.impurity)
            print('left/right', self.criterion.impurity_left,
                   self.criterion.impurity_right)
            # check the impurity improved
            if current_impurity < self.split_record.impurity:
                # Update all the record for the current record
                self.split_record.threshold = (
                    (self.X[sample_idx, feat_i] +
                     self.X[self.prev_idx, feat_i]) / 2.0)
                self.split_record.impurity = current_impurity
                self.split_record.pos = self.prev_idx
                # Update the statistic of the future children
                self.split_record.sum_left_residual = self.criterion.sum_left
                self.split_record.sq_sum_left_residual = self.criterion.sq_sum_left
                self.split_record.n_left_samples = self.criterion.n_left_samples
                self.split_record.weighted_left_samples = self.criterion.weighted_n_left
                # right node
                self.split_record.sum_right_residual = self.criterion.sum_right
                self.split_record.sq_sum_right_residual = self.criterion.sq_sum_right
                self.split_record.n_right_samples = self.criterion.n_right_samples
                self.split_record.weighted_right_samples = self.criterion.weighted_n_right

        self.prev_idx = sample_idx
