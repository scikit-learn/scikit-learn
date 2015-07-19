# Title          : discretization.py
# Author         : Henry Lin <hlin117@gmail.com>
# License        : BSD 3 clause
#==============================================================================

from __future__ import division
from math import log
import numpy as np
from scipy.stats import entropy
from ..base import BaseEstimator
from ..utils import check_array, check_consistent_length, column_or_1d
import abc

class KBins(TransformerMixin):
    """TODO
    """

    def __init__(self, numBins, continuous_features=None):
        if numBins < 2:
            raise ValueError("KBins discretizer received an invalid number "
                             "of bins")
        self.numBins = numBins
        self.continuous_features_ = continuous_features

    def fit(self, X, y=None):
        """Finds the intervals of interest from the input data.

        Parameters
        ----------
        X : The array containing features to be discretized. Continuous
        features should be specified by the `continuous_features`
        attribute.

        y : A list or array of class labels corresponding to `X`.
        """
        X = check_array(X, force_all_finite=True, ensure_2d=False)
        if self.continuous_features_ is None:
            self.continuous_features_ = range(X.shape[1])

        self.cut_points_ = dict()

        for index, col in enumerate(X.T):
            if index not in self.continuous_features_:
                continue
            cut_points = self._split_intervals(col, y)
            self.cut_points_[index] = cut_points
        return self

    def transform(self, X, y=None):
        """Converts the continuous features in X into integers from
        0... k-1 (`k` is the number of intervals the discretizer created
        from a given continuous feature.)
        """
        if self.continuous_features_ is None:
            raise ValueError("You must fit the discretizer before "
                             "transforming data.")

        output = X.copy()
        for i in self.continuous_features_:
            output[:, i] = np.searchsorted(self.cut_points_[i], X[:, i])
        return output

    def _split_intervals(self, col, y=None):
        """Fixed length binning.

        Returns
        -------
        A python list of cut points, representing the boundaries
        of each interval.
        """
        # Create the intervals
        maxVal = max(col)
        minVal = min(col)
        spacing = (maxVal - minVal) / (self.numBins)
        cut_points = np.array([(minVal + spacing * i)
                                for i in range(1, self.numBins)])
        return cut_points
