# Title          : discretization.py
# Author         : Henry Lin <hlin117@gmail.com>
# License        : BSD 3 clause
#==============================================================================

from __future__ import division
from math import log
import numpy as np
from scipy.stats import entropy
from sklearn.base import BaseEstimator

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
        found to be zero before `min_depth`, the algorithm will automatically
        stop.

        Defaults to 0.

    num_classes_ : The number of classes the discretizer was fit to in `y`.

    intervals_ : A dictionary from indices to lists. Each list contains
        tuples of the form (interval, cat), where `interval` is a pair
        `(a, b)`, and `cat` is an integer from 0... k-1 (`k` is the number
        of intervals a continuous attribute has.)

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
    continuous interval (a, b], where `a` can be `float("-inf")`, and `b` can
    be `float("inf")`.
    """

    def __init__(self, continuous_features=None, min_depth=0):
        self.continuous_features_ = continuous_features
        self.min_depth = min_depth

    def fit(self, X, y):
        """Finds the intervals of interest from the input data.
        """
        if type(X) is list:
            X = np.array(X)
        if type(y) is list:
            y = np.array(y)
        if self.continuous_features_ is None:
            self.continuous_features_ = range(X.shape[1])
        if(len(X.shape) != 2):
            raise ValueError("MDLP can ony be applied to input ndarrays of "
                             "size 2.")
        if(X.shape[0] != y.shape[0]):
            raise ValueError("Number of samples in X does not match number "
                             "of targets in y")

        self.num_classes_ = set(y)
        self.intervals_ = dict()

        for index, col in enumerate(X.T):
            if index not in self.continuous_features_:
                continue
            intervals = self._mdlp(col, y)
            intervals.sort(key=lambda interval: interval[0])
            self.intervals_[index] = [(interval, i)
                                     for i, interval in enumerate(intervals)]
        return self

    def transform(self, X, y=None):
        """Converts the continuous features in X into integers from
        0... k-1 (`k` is the number of intervals the discretizer created
        from a given continuous feature.)
        """
        assert self.continuous_features_ is not None, "You must fit the object " \
                                                    "before transforming data."
        discretized = list()
        for index, col in enumerate(X.T):
            if index not in self.continuous_features_:
                discretized.append(col.reshape(len(X), 1))
            transformed = self.cts2cat(index, col)
            discretized.append(transformed.reshape(len(X), 1))
        return np.hstack(discretized)

    def fit_transform(self, X, Y):
        """Performs the fitting and the data transformations.
        """
        self.fit(X, Y)
        return self.transform(X)

    def cat2intervals(self, X, index):
        """Converts a categorical feature into a list of intervals.
        This is the naive linear search method, but in practice, because
        there will not  be many intervals to search through (usually at
        most 5), this works effectively without an interval tree.
        """
        if index not in self.continuous_features_:
            raise ValueError("The input index {0} is not listed as a "
                             "continuous feature".format(index))
        col = X.T[index]
        mappings = self.intervals_[index]
        out_intervals = list()
        for cat in col:
            converted = False
            for interval, currcat in mappings:
                if cat == currcat:
                    out_intervals.append(interval)
                    converted = True
                    break
            if not converted:
                raise ValueError("The categorical value {0} was not found "
                                 "in this column.".format(cat))
        return out_intervals

    def bin(self, index, attr):
        """Bins an attribute into its appropriate interval.
        """
        mappings = self.intervals_[index]
        for (a, b), _ in mappings:
            if a < attr <= b: return a, b
        raise ValueError("Numeric value did not fit in any interval")

    def cts2cat(self, index, col):
        """Converts each continuous feature into a categorical feature.
        """
        mappings = self.intervals_[index]
        categorical = list()
        for attr in col:
            for (a, b), cat in mappings:
                if a < attr <= b:
                    categorical.append(cat)
                    break
        assert len(categorical) == len(col)
        return np.array(categorical)

    def _mdlp(self, x, y):
        """
        *attributes*: A numpy 1 dimensional ndarray

        *y*: A python list of numeric class labels.

        Returns a list of intervals, based upon the MDLP stopping criterion.
        """

        order = np.argsort(x)
        x = x[order]
        y = y[order]

        intervals = list()

        def get_cut(ind):
            return (x[ind-1] + x[ind]) / 2

        # Now we do a depth first search to create intervals
        num_samples = len(x)
        search_intervals = list()
        search_intervals.append((0, num_samples, 0))
        while len(search_intervals) > 0:
            start, end, depth = search_intervals.pop()

            k = self._find_cut(y, start, end)

            # Need to see whether the "front" and "back" of the intervals need
            # to be float("-inf") or float("inf")
            if (k == -1) or (self._reject_split(y, start, end, k) and
                             depth >= self.min_depth):
                front = float("-inf") if (start == 0) else get_cut(start)
                back = float("inf") if (end == num_samples) else get_cut(end)

                # Corner case: If front == back, don't add it in
                if front < back:
                    intervals.append((front, back))
                continue

            search_intervals.append((start, k, depth + 1))
            search_intervals.append((k, end, depth + 1))

        return intervals

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

        def partition_entropy(ind):
            first_half = (ind - start) / length * \
                    MDLP._slice_entropy(y, start, ind)[0]
            second_half = (end - ind) / length * \
                    MDLP._slice_entropy(y, ind, end)[0]
            return first_half + second_half

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

            curr_entropy = partition_entropy(ind)
            if prev_entropy > curr_entropy:
                prev_entropy = curr_entropy
                k = ind

        return k  # NOTE: k == -1 if there is no good cut
