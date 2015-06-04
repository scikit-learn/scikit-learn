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

class MDLP(BaseEstimator):
    """Implements the MDLP discretization criterion from Usama Fayyad's
    paper "Multi-Interval Discretization of Continuous-Valued Attributes
    for Classification Learning". Given the class labels for each sample,
    this transformer attempts to discretize a continuous attribute by
    minimizing the entropy at each interval.

    Attributes
    ----------
    min_depth : The minimum depth of the interval splitting. Overrides the
        MDLP stopping criterion. If the entropy at a given interval is
        found to be zero before `min_depth`, the algorithm will stop.

        Defaults to 0.

    num_classes_ : The number of classes the discretizer was fit to in `y`.

    cut_points_ : A dictionary mapping indices to a numpy array. Each
    numpy array is a sorted list of cut points found from discretization.

    continuous_features_ : A list of indices indicating which columns should
        be discretized.

        Defaults to `range(X.shape[1])` if not specified.

    Examples
    --------

    >>> from sklearn.preprocessing import MDLP
    >>> from sklearn.datasets import load_iris
    >>> iris = load_iris()
    >>> X = iris.data
    >>> y = iris.target
    >>> mdlp = MDLP()
    >>> conv_X = mdlp.fit_transform(X, y)

    `conv_X` will be the same shape as `X`, except it will contain
    integers instead of continuous attributes representing the results of
    the discretization process.

    To retrieve the explicit intervals of the discretization of, say, the
    third column (index 2), one can do

    >>> mdlp.cat2intervals(conv_X, 2)

    which would return a list of tuples `(a, b)`. Each tuple represents the
    contnuous interval (a, b], where `a` can be `float("-inf")`, and `b` can
    be `float("inf")`.
    """

    def __init__(self, continuous_features=None, min_depth=0):
        self.continuous_features_ = continuous_features
        self.min_depth = min_depth

    def fit(self, X, y):
        """Finds the intervals of interest from the input data.

        Parameters
        ----------
        X : The array containing features to be discretized. Continuous
            features should be specified by the `continuous_features`
            attribute.

        y : A list or array of class labels corresponding to `X`.
        """
        X = check_array(X, force_all_finite=True, estimator="MDLP discretizer")
        y = column_or_1d(y)
        check_consistent_length(X, y)
        if self.continuous_features_ is None:
            self.continuous_features_ = range(X.shape[1])

        self.num_classes_ = set(y)
        self.cut_points_ = dict()

        for index, col in enumerate(X.T):
            if index not in self.continuous_features_:
                continue
            cut_points = self._mdlp(col, y)
            self.cut_points_[index] = np.array(cut_points)
        return self

    def transform(self, X, y=None):
        """Converts the continuous features in X into integers from
        0... k-1 (`k` is the number of intervals the discretizer created
        from a given continuous feature.)
        """
        if self.continuous_features_ is None:
            raise ValueError("You must fit the MDLP discretizer before "
                             "transforming data.")

        output = X.copy()
        for i in self.continuous_features_:
            output[:, i] = np.searchsorted(self.cut_points_[i], X[:, i])
        return output

    def fit_transform(self, X, Y):
        """Performs the fitting and the data transformations.
        """
        self.fit(X, Y)
        return self.transform(X)

    def cat2intervals(self, X, index):
        """Converts a categorical feature into a list of intervals.
        """
        cp_indices = X.T[index]
        return self._assign_intervals(cp_indices, index)

    def _assign_intervals(self, cp_indices, index):
        """Assigns the cut point indices `cp_indices` (representing
        categorical features) into a list of intervals.
        """

        non_zero_mask = cp_indices[cp_indices - 1 != -1].astype(int) - 1
        fronts = np.zeros(cp_indices.shape)
        fronts[cp_indices == 0] = float("-inf")
        fronts[cp_indices != 0] = self.cut_points_[index][non_zero_mask]

        numCuts = len(self.cut_points_[index])
        non_numCuts_mask = cp_indices[cp_indices != 2].astype(int)
        backs = np.zeros(cp_indices.shape)
        backs[cp_indices == numCuts] = float("inf")
        backs[cp_indices != numCuts] = self.cut_points_[index][non_numCuts_mask]

        return [(front, back) for front, back in zip(fronts, backs)]


    def cts2intervals(self, col, index):
        """Bins a continuous feature into its appropriate interval.

        TODO: There's probably a more elegant way of going about this
        when `col` is a scalar.
        """
        cp_indices = self.cts2cat(col, index)
        if cp_indices is not np.array:
            cp_indices = np.array([cp_indices])
        return self._assign_intervals(cp_indices, index)

    def cts2cat(self, col, index):
        """Converts each continuous feature from index `index` into
        a categorical feature from the input column `col`.
        """
        return np.searchsorted(self.cut_points_[index], col)

    def _mdlp(self, col, y):
        order = np.argsort(col)
        col = col[order]
        y = y[order]

        cut_points = set()

        def get_cut(ind):
            return (col[ind-1] + col[ind]) / 2

        # Now we do a depth first search to create cut_points
        num_samples = len(col)
        search_intervals = list()
        search_intervals.append((0, num_samples, 0))
        while len(search_intervals) > 0:
            start, end, depth = search_intervals.pop()

            k = self._find_cut(y, start, end)

            # Need to see whether the "front" and "back" of the interval need
            # to be float("-inf") or float("inf")
            if (k == -1) or (depth >= self.min_depth and
                             self._reject_split(y, start, end, k)):
                front = float("-inf") if (start == 0) else get_cut(start)
                back = float("inf") if (end == num_samples) else get_cut(end)

                if front == back: continue  # Corner case
                if front != float("-inf"): cut_points.add(front)
                if back != float("inf"): cut_points.add(back)
                continue

            search_intervals.append((start, k, depth + 1))
            search_intervals.append((k, end, depth + 1))

        cut_points = list(cut_points)
        cut_points.sort()
        return cut_points

    @staticmethod
    def _slice_entropy(y, start, end):
        """Returns the entropy of the given slice of y. Also returns the
        number of classes within the interval.
        """
        counts = np.bincount(y[start:end])
        vals = counts / (end - start)
        return entropy(vals), np.sum(vals != 0)

    @staticmethod
    def _reject_split(y, start, end, k):
        """Using the minimum description length principal, determines
        whether it is appropriate to stop cutting.
        """

        N = end - start
        entropy1, k1 = MDLP._slice_entropy(y, start, k)
        entropy2, k2 = MDLP._slice_entropy(y, k, end)
        whole_entropy, k = MDLP._slice_entropy(y, start, end)

        # Calculate the final values
        gain = whole_entropy - 1 / N * ((start - k) * entropy1 +
                                        (end - k) * entropy2)
        entropy_diff = k * whole_entropy - k1 * entropy1 - k2 * entropy2
        delta = log(3**k - 2) - entropy_diff

        return gain <= 1 / N * (log(N - 1) + delta)

    @staticmethod
    def _find_cut(y, start, end):
        """Finds the best cut between the specified interval.

        Returns a candidate cut point index k. k == -1 if there is no good cut,
        otherwise k is in {1, ..., len(x) - 2}
        """

        # Want to probe for the best partition _entropy in a "smart" way
        # Input is the splitting index, to create partitions [start, ind)
        # and [ind, end).
        length = end - start

        """
        For each iteration, we'll have start != ind and ind != end
        We'll also have the length of both partitions at least 1

        TODO: The smarter method would be to do a one pass to collect
        the number of positives and negatives for each cut point, and
        perform entropy calculations that way. This method takes n^2 time per
        iteration because of `_slice_entropy()`.
        """
        prev_entropy = float("inf")
        k = -1
        for ind in range(start + 1, end - 1):

            # I choose not to use a `min` function here for this optimization.
            if y[ind-1] == y[ind]:
                continue

            # Finds the partition entropy, and see if this entropy is minimum
            first_half = (ind - start) / length * \
                    MDLP._slice_entropy(y, start, ind)[0]
            second_half = (end - ind) / length * \
                    MDLP._slice_entropy(y, ind, end)[0]
            curr_entropy = first_half + second_half

            if prev_entropy > curr_entropy:
                prev_entropy = curr_entropy
                k = ind

        return k  # NOTE: k == -1 if there is no good cut
