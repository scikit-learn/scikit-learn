#title          :discretization.py
#description    :Implementation of a few discretization algorithms
#author         :Henry Lin <hlin117@gmail.com>
#==============================================================================

from __future__ import division
from math import log
import numpy as np
from scipy.stats import entropy


def log2(x):
    return log(x, 2) if x > 0 else 0


class MDLP(object):
    """Implements the MDLP discretization criterion from Usama Fayyad's
    paper "Multi-Interval Discretization of Continuous-Valued Attributes
    for Classification Learning"
    """

    def __init__(self, **params):
        self.intervals_ = dict()
        self.num_classes_ = None
        self.set_params(**params)

    def set_params(self, **params):
        self.continuous_columns = params.get("continous_columns")
        self.min_depth = params.get("min_depth")

    def fit(self, X, Y):
        """Finds the intervals of interest from the input data.
        """
        if type(X) is list:
            X = np.array(X)
        if self.continuous_columns is None:
            self.continuous_columns = range(X.shape[1])
        assert len(X.shape) == 2, "MDLP can ony be applied to ndarrays of " \
                                  "size 2."
        self.num_classes_ = set(Y)

        for index, col in enumerate(X.T):
            if index not in self.continuous_columns:
                continue
            intervals = self._mdlp(col, Y, self.min_depth)
            intervals.sort(key=lambda interval: interval[0])
            self.intervals_[index] = [(interval, i)
                                     for i, interval in enumerate(intervals)]

    def transform(self, X):
        """Converts the continuous values in X into ascii character
        values. The mapping is defined by self.intervals.
        """
        assert self.continuous_columns is not None, "You must fit the object " \
                                                    "before transforming data."
        discretized = list()
        for index, col in enumerate(X.T):
            if index not in self.continuous_columns:
                discretized.append(col.reshape(len(X), 1))
            transformed = self.cts2cat(index, col)
            discretized.append(transformed.reshape(len(X), 1))
        return np.hstack(discretized)

    def fit_transform(self, X, Y):
        """Performs the fitting and the data transformations.
        """
        self.fit(X, Y)
        return self.transform(X)

    def cat2intervals(self, index, col):
        """Converts a categorical column into an interval. This is the naive
        linear search method, but in practice, because there won't be many
        intervals to search through (usually at most 5), this works effectively
        without an interval tree.
        """
        mappings = self.intervals_[index]
        out_intervals = list()
        for cat in col:
            for interval, currcat in mappings:
                if cat == currcat:
                    out_intervals.append(interval)
                    break
        return out_intervals

    def bin(self, index, attr):
        """Bins an attribute into its appropriate interval.
        """
        mappings = self.intervals_[index]
        for (a, b), _ in mappings:
            if a < attr <= b: return a, b
        raise ValueError("Numeric value did not fit in any interval")

    def cts2cat(self, index, col):
        """Converts each continuous value into a categorical value.
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

    def _mdlp(self, attributes, Y, min_depth):
        """
        *attributes*: A numpy 1 dimensional ndarray

        *Y*: A python list of numeric class labels.

        Returns a list of intervals, based upon the MDLP stopping criterion.
        """

        order = np.argsort(attributes)
        x = attributes[order]
        y = Y[order]

        intervals = list()

        def get_cut(ind):
            return ((x[ind-1] + x[ind]) / 2)

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
                             depth >= min_depth):
                front = float("-inf") if (start == 0) else get_cut(start)
                back = float("inf") if (end == num_samples) else get_cut(end)

                # Corner case: If front == back, don't add it in
                if front < back:
                    intervals.append((front, back))
                continue

            search_intervals.append((start, k, depth + 1))
            search_intervals.append((k, end, depth + 1))

        return intervals

    def _slice_entropy(self, y, start, end):
        """Returns the entropy of the given slice of y. Also returns the
        number of classes within the interval.
        """
        counts = np.bincount(y[start:end])
        vals = counts / (end - start)
        return entropy(vals, base=2), np.sum(vals != 0)

    def _reject_split(self, y, start, end, k):
        """Using the minimum description length principal, determines
        whether it is appropriate to stop cutting.
        """

        N = end - start
        entropy1, k1 = self._slice_entropy(y, start, k)
        entropy2, k2 = self._slice_entropy(y, k, end)
        whole_entropy, k = self._slice_entropy(y, start, end)

        # Calculate the final values
        gain = whole_entropy - 1 / N * ((start - k) * entropy1 +
                                        (end - k) * entropy2)
        entropy_diff = k * whole_entropy - k1 * entropy1 - k2 * entropy2
        delta = log2(3**k - 2) - entropy_diff

        return gain <= 1 / N * (log2(N - 1) + delta)

    def _find_cut(self, y, start, end):
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
                    self._slice_entropy(y, start, ind)[0]
            second_half = (end - ind) / length * \
                    self._slice_entropy(y, ind, end)[0]
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
        for ind in xrange(start + 1, end - 1):

            # I choose not to use a `min` function here for this optimization.
            if y[ind-1] == y[ind]:
                continue

            curr_entropy = partition_entropy(ind)
            if prev_entropy > curr_entropy:
                prev_entropy = curr_entropy
                k = ind

        return k  # NOTE: k == -1 if there is no good cut
