from __future__ import division, print_function

from .criterion import NewCriterion
from .criterion import NewSplitRecord

FEATURE_THRESHOLD = 1e-7

class NewSplitter(object):
    """New type of splitter driven by data
    """

    def __init__(self, X, y, sample_weight, weighted_n_total_samples, nid,
                 p_split_record, min_samples_leaf, min_weight_leaf):
        self.X = X
        self.y = y
        self.sample_weight = sample_weight
        self.weighted_n_total_samples = weighted_n_total_samples
        self.nid = nid

        # get the parent criterion
        self.p_split_record = p_split_record

        # parameters for early stop of split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_leaf = min_weight_leaf

        # split record for the current node
        self.split_record = NewSplitRecord()
        self.split_record.nid = nid

        # create the criterion for the current splitter
        self.criterion = NewCriterion(self.X, self.y, self.sample_weight,
                                      self.weighted_n_total_samples,
                                      self.p_split_record)

    def reset(self, X):
        """Reset a splitter with a new samples set.
        It could correspond to a new feature to be scanned.
        """
        self.X = X

        # reinitialize the split record
        self.split_record.__init__()
        self.split_record.nid = self.nid

        # reinitialize the criterion
        self.criterion.__init__(self.X, self.y, self.sample_weight,
                                self.weighted_n_total_samples,
                                self.p_split_record)

    @property
    def split_record(self):
        """Get the split_record for that node"""
        return self.split_record

    def node_evaluate_split(self, sample_idx):
        """Update the impurity and check the corresponding split should be
        kept.
        """
        # check that the sample value are different enough
        if self.X[sample_idx] <= self.X[self.split_record.position]:

            # make an update of the statistics
            self.criterion.update_stats(sample_idx)

            # check that there is enough samples to make the proper split
            if (self.criterion.n_left_samples < self.min_samples_leaf or
                    self.criterion.n_right_samples < self.min_samples_leaf):
                # early stopping
                return

            # check that the weight are enough important to make a proper split
            if (self.criterion.weighted_n_left < self.min_weight_leaf or
                    self.criterion.weighted_n_right < self.min_weight_leaf):
                # early stopping
                return

            # compute the impurity improvement
            current_impurity = self.criterion.impurity_improvement()

            # check the impurity improved
            if current_impurity > self.split_record.impurity:
                self.split_record.threshold = (
                    (self.X[sample_idx] +
                     self.X[self.split_record.position]) / 2.0)
                self.split_record.impurity = current_impurity
                self.split_record.position = sample_idx
