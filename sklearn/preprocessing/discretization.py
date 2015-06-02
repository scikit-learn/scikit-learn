#title          :discretization.py
#description    :Implementation of a few discretization algorithms
#author         :Henry Lin <hlin117@gmail.com>
#version        :0.0.1
#python_version :2.7.6
#==============================================================================

from collections import Counter
from math import log
import numpy as np
from itertools import izip, islice


def log2(x):
    return log(x, 2) if x > 0 else 0


class MDLP(object):
    """Implements the MDLP algorithm from Usama Fayyad's paper
    "Multi-Interval Discretization of Continuous-Valued Attributes for
    Classification Learning"
    """

    def __init__(self):
        self.intervals = dict()
        self.continuous_columns = None

    def fit(self, X, Y, continuous_columns=[], min_depth=None):
        """Finds the intervals of interest from the input data.
        """
        self.continuous_columns = continuous_columns
        if type(X) is list:
            X = np.array(X)
        if len(self.continuous_columns) == 0:
            self.continuous_columns = range(X.shape[1])
        assert len(X.shape) == 2, "MDLP can ony be applied to ndarrays of " \
                                  "size 2."
        for index, col in enumerate(X.T):
            if index not in self.continuous_columns:
                continue
            intervals = self._mdlp(col, Y, min_depth)
            intervals.sort(key=lambda (a, b): a)
            self.intervals[index] = [(interval, i)
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

    def fit_transform(self, X, Y, continuous_columns=[], min_depth=None):
        """Performs the fitting and the data transformations.
        """
        self.fit(X, Y, continuous_columns, min_depth)
        return self.transform(X)

    def cat2intervals(self, index, col):
        """Converts a categorical column into an interval. This is the naive
        linear search method, but in practice, because there won't be many
        intervals to search through (usually at most 5), this works effectively
        without an interval tree.
        """
        mappings = self.intervals[index]
        out_intervals = list()
        for cat in col:
            for interval, currcat in mappings:
                if cat == currcat:
                    out_intervals.append(interval)
        return out_intervals

    def bin(self, index, attr):
        """Bins an attribute into its appropriate interval.
        """
        mappings = self.intervals[index]
        for (a, b), _ in mappings:
            if a < attr <= b: return (a, b)
        raise ValueError("Numeric value did not fit in any interval")

    def cts2cat(self, index, col):
        """Converts each continuous value into a categorical value.
        """
        mappings = self.intervals[index]
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
        *attributes*: A python list of continuous values.

        *Y*: A python list of 0's and 1's

        Returns a list of tuples representing intervals.
        """

        # attr_list contains triples, (attribute, y class)
        attr_list = list([attr, y] for attr, y in izip(attributes, Y))
        attr_list.sort(key=lambda (attr, y): attr)
        intervals = list()

        def get_cut(ind):
            return ((attr_list[ind-1][0] + attr_list[ind][0]) / 2.0)

        # Now we do a depth first search to create intervals
        # attr_list will be changed, so each attr_list[0] will contain a tuple
        # rather than a continuous value.
        search_intervals = list()
        search_intervals.append((0, len(attr_list), 0))
        while len(search_intervals) > 0:
            start, end, depth = search_intervals.pop()

            k = self._find_cut(attr_list, start, end)

            # Need to see whether the "front" and "back" of the intervals need
            # to be infinity intervals
            if (k == -1) or (self._reject_split(attr_list, start, end, k) and
                             depth >= min_depth):
                front = float("-inf") if (start == 0) else get_cut(start)
                back = float("inf") if (end == len(attr_list)) else get_cut(end)

                # Corner case: If front == back, don't add it in
                if front < back:
                    intervals.append((front, back))
                continue

            search_intervals.append((start, k, depth + 1))
            search_intervals.append((k, end, depth + 1))

        # Place items back in the order they appeared, and return
        return intervals

    def _reject_split(self, attr_list, start, end, k):
        """Using the minimum description length principal, determines
        whether it is appropriate to stop cutting.

        *attr_list*: A list of sorted continuous attributes that we would
        like to find a good cut point of. Each item is a triple:
        (cts_val, y), where y is a class label.

        "pos" here reflects large error. "neg" here reflect small error.
        """

        cat_counter1 = Counter()
        cat_counter2 = Counter()

        for (attr, y) in islice(attr_list, start, k):
            cat_counter1[y] += 1
        for (attr, y) in islice(attr_list, k, end):
            cat_counter2[y] += 1

        big_counter = {y: (cat_counter1[y] + cat_counter2[y])
                       for y in set(cat_counter1.keys() + cat_counter2.keys())}

        entropy1 = self._counter_entropy(cat_counter1, start, k)
        entropy2 = self._counter_entropy(cat_counter2, k, end)
        big_entropy = self._counter_entropy(big_counter, start, end)

        N = end - start

        # Calculate the final values
        entropy_diff = len(big_counter) * big_entropy - \
                len(cat_counter1) * entropy1 - len(cat_counter2) * entropy2
        delta = log2(3**len(big_counter) - 2) - entropy_diff
        gain = big_entropy - 1 / float(N) * ((k - start) * entropy1 +
                                             (end - k) * entropy2)

        return gain <= 1 / float(N) * (log2(N - 1) + delta)

    def _find_cut(self, attr_list, start, end):
        """Finds the best cut between the specified interval.

        *attr_list*: A list of sorted continuous attributes that we would
        like to find a good cut point of.

        Each item is a triple: (cts_val, y, i),
        where y is a class label, and i is the index in the original dataset.

        *start*: The beginning of the list, inclusive.

        *end*: The end of the list, _exclusive_

        Returns a continuous value cut_point, and an index k such that
        cut_point < attr_list[k][0]. k == -1 if there is no good cut,
        otherwise k is in {1, ..., len(attr_list) - 2}
        """

        # Want to probe for the best partition _entropy in a "smart" way
        # Input is the splitting index, to create partitions [start, ind)
        # and [ind, end).
        length = end - start

        def partition_entropy(ind):
            first_half = (ind - start) / float(length) * \
                    self._entropy(attr_list, start, ind)
            second_half = (end - ind) / float(length) * \
                    self._entropy(attr_list, ind, end)
            return first_half + second_half

        """
        For each iteration, we'll have start != ind and ind != end
        We'll also have the length of both partitions at least 1

        TODO: The smarter method would be to do a one pass to collect
        the number of positives and negatives for each point, and
        perform _entropy calculations that way. This method takes n^2 time per
        iteration.
        """
        prev_entropy = float("inf")
        k = -1
        for ind in xrange(start + 1, end - 1):

            # I choose not to use a `min` function here for this optimization.
            if attr_list[ind-1][1] == attr_list[ind][1]:
                continue

            curr_entropy = partition_entropy(ind)
            if prev_entropy > curr_entropy:
                prev_entropy = curr_entropy
                k = ind

        return k  # NOTE: k == -1 if there is no good cut

    def _entropy(self, attr_list, start, end):
        """Calculates the _entropy of the subset contained in the interval
        start and end (exclusive)
        """
        cat_counter = Counter()
        for (attr, y) in islice(attr_list, start, end):
            cat_counter[y] += 1

        return float(self._counter_entropy(cat_counter, start, end))

    def _counter_entropy(self, counter, start, end):
        """A helper function, specifically if you have counted already.
        """
        length = end - start
        return -sum(count / float(length) * log2(count / float(length))
                    for count in counter.itervalues())
