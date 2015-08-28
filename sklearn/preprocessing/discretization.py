# Title          : discretization.py
# Author         : Henry Lin <hlin117@gmail.com>
# License        : BSD 3 clause
#==============================================================================

from __future__ import division
import numpy as np
from ..base import BaseEstimator
from sklearn.base import TransformerMixin
from ..utils import check_array, check_X_y, column_or_1d
from _mdlp import _find_cut, _reject_split

class MDLP(BaseEstimator, TransformerMixin):
    """Implements the MDLP discretization criterion from Usama Fayyad's
    paper "Multi-Interval Discretization of Continuous-Valued
    Attributes for Classification Learning". Given the class labels
    for each sample, this transformer attempts to discretize a
    continuous attribute by minimizing the entropy at each interval.

    Attributes
    ----------
    min_depth : The minimum depth of the interval splitting. Overrides
        the MDLP stopping criterion. If the entropy at a given interval
        is found to be zero before `min_depth`, the algorithm will stop.

        Defaults to 0.

    cut_points_ : A dictionary mapping indices to a numpy array. Each
        numpy array is a sorted list of cut points found from
        discretization.

    dimensions_ : Describes whether `X` is a 2-D or 1-D array.

    continuous_features_ : A list of indices indicating which columns
        should be discretized.

        If `X` is a 1-D array, then should be None. Otherwise, for a
        2-D array, defaults to `range(X.shape[1])`.

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
    integers instead of continuous attributes representing the results
    of the discretization process.

    To retrieve the explicit intervals of the discretization of, say,
    the third column (index 2), one can do

    >>> mdlp.cat2intervals(conv_X, 2)

    which would return a list of tuples `(a, b)`. Each tuple represents
    the contnuous interval (a, b], where `a` can be `float("-inf")`,
    and `b` can be `float("inf")`.
    """

    def __init__(self, continuous_features=None, min_depth=0, shuffle=True):
        self.continuous_features_ = continuous_features
        self.min_depth = min_depth
        self.shuffle = shuffle
        self.cut_points_ = None
        self.dimensions_ = None

    def fit(self, X, y):
        """Finds the intervals of interest from the input data.

        Parameters
        ----------
        X : The array containing features to be discretized. Continuous
            features should be specified by the `continuous_features`
            attribute if `X` is a 2-D array.

        y : A list or array of class labels corresponding to `X`.
        """
        self.dimensions_ = len(X.shape)

        if self.dimensions_ > 2:
            raise ValueError("Invalid input dimension for `X`. Input shape is"
                             "{0}".format(X.shape))

        X = check_array(X, force_all_finite=True, estimator="MDLP discretizer",
                        ensure_2d=False)
        y = column_or_1d(y)
        y = check_array(y, ensure_2d=False, dtype=int)
        X, y = check_X_y(X, y)

        if self.dimensions_ == 2:
            if self.continuous_features_ is None:
                self.continuous_features_ = range(X.shape[1])

            self.cut_points_ = dict()

            for index, col in enumerate(X.T):
                if index not in self.continuous_features_:
                    continue
                cut_points = self._mdlp(col, y)
                self.cut_points_[index] = np.array(cut_points)
        else:
            # TODO: Throw warning if self.continuous_features_ is not None
            self.continuous_features_ = None
            cut_points = self._mdlp(X, y)
            self.cut_points_ = np.array(cut_points)

        return self

    def transform(self, X, y=None):
        """Converts the continuous features in X into integers from
        0... k-1 (`k` is the number of intervals the discretizer created
        from a given continuous feature.)
        """
        if self.cut_points_ is None:
            raise ValueError("You must fit the MDLP discretizer before "
                             "transforming data.")
        if self.dimensions_ == 1:
            output = np.searchsorted(self.cut_points_, X)
        else:
            output = X.copy()
            for i in self.continuous_features_:
                output[:, i] = np.searchsorted(self.cut_points_[i], X[:, i])
        return output

    def cat2intervals(self, X, index=None):
        """Converts a categorical feature into a list of intervals.
        """
        # TODO: Throw warning if `self.dimensions_` == 1 and index is not None
        if self.dimensions_ == 1:
            return self._assign_intervals(X, index)
        elif self.dimensions_ == 2 and index is None:
            raise ValueError("Index of `X` to be discretized needs to be "
                             "specified.")
        else:
            cp_indices = X.T[index]
            return self._assign_intervals(cp_indices, index)

    def _assign_intervals(self, cp_indices, index):
        """Assigns the cut point indices `cp_indices` (representing
        categorical features) into a list of intervals.
        """

        # Case for a 1-D array
        if self.dimensions_ == 1:
            cut_points = self.cut_points_
        else:
            cut_points = self.cut_points_[index]

        non_zero_mask = cp_indices[cp_indices - 1 != -1].astype(int) - 1
        fronts = np.zeros(cp_indices.shape)
        fronts[cp_indices == 0] = float("-inf")
        fronts[cp_indices != 0] = cut_points[non_zero_mask]

        numCuts = len(cut_points)
        backs = np.zeros(cp_indices.shape)
        non_numCuts_mask = cp_indices[cp_indices != numCuts].astype(int)
        backs[cp_indices == numCuts] = float("inf")
        backs[cp_indices != numCuts] = cut_points[non_numCuts_mask]

        return [(front, back) for front, back in zip(fronts, backs)]

    def cts2cat(self, col, index=None):
        """Converts each continuous feature from index `index` into
        a categorical feature from the input column `col`.
        """
        if self.dimensions_ == 1:
            return np.searchsorted(self.cut_points_, col)
        if self.dimensions_ == 2 and index is None:
            raise ValueError("Index of `X` to be discretized needs to be "
                             "specified.")
        return np.searchsorted(self.cut_points_[index], col)

    def cts2intervals(self, col, index=None):
        """Bins a continuous feature into its appropriate interval.
        """
        cp_indices = self.cts2cat(col, index)

        # TODO: Below occurs if `col` is a number. There's probably a more
        # elegant way of handling this case.
        if cp_indices is not np.array:
            cp_indices = np.array([cp_indices])
        return self._assign_intervals(cp_indices, index)

    def _mdlp(self, col, y):

        # Shuffle array, and then reorder them
        if self.shuffle:
            shuffled_order = np.random.permutation(len(y))
            col = col[shuffled_order]
            y = y[shuffled_order]

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

            k = _find_cut(y, start, end)

            # Need to see whether the "front" and "back" of the interval need
            # to be float("-inf") or float("inf")
            if (k == -1) or (depth >= self.min_depth and
                             _reject_split(y, start, end, k)):
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

