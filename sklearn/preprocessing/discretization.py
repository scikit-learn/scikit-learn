# Title          : discretization.py
# Author         : Henry Lin <hlin117@gmail.com>
# License        : BSD 3 clause
#==============================================================================

from __future__ import division
import numpy as np
from ..externals import six
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import column_or_1d
from sklearn.utils.validation import check_is_fitted

xrange = six.moves.xrange

__all__ = [
    "FixedWidthDiscretizer"
]


class FixedWidthDiscretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    n_bins : int (default=2)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.

    Attributes
    ----------
    min_ : float
        The minimum value of the input data.

    max_ : float
        The maximum value of the input data.

    cut_points_ : array, shape [numBins - 1]
        Contains the boundaries for which the data lies. Each interval
        has an open left boundary, and a closed right boundary.

    spacing_ : float
        The spacing between each one of the cut points. Calculated as

            (max - min) / n_bins.

    Example
    -------
    >>> from sklearn.preprocessing import FixedWidthDiscretizer
    >>> X = [0, 1, 2, 3, 4, 5, 6]
    >>> discretizer = FixedWidthDiscretizer(n_bins=3)
    >>> discretizer.fit(X)
    FixedWidthDiscretizer(n_bins=3)
    >>> discretizer.cut_points_
    array([2., 4.])
    >>> discretizer.transform(X)
    array([0, 0, 0, 1, 1, 2, 2])
    """

    def __init__(self, n_bins=2):
        if n_bins < 2:
            raise ValueError("FixedWidthDiscretizer received an invalid number "
                             "of bins")
        self.n_bins = n_bins

        # Attributes
        self.min_ = None
        self.max_ = None
        self.cut_points_ = None
        self.spacing_ = None

    def fit(self, X, y=None):
        """Finds the intervals of interest from the input data.

        Parameters
        ----------
        X : array-like, shape [n_samples]
            The array containing continuous features to be discretized.
            Input must be 1d arrays.

        """

        X = column_or_1d(X)

        self.min_ = X.min()
        self.max_ = X.max()
        self.spacing_ = (self.max_ - self.min_) / self.n_bins
        cut_points = (self.min_ + self.spacing_ * i
                      for i in xrange(1, self.n_bins))
        self.cut_points_ = np.array(list(cut_points))
        return self

    def transform(self, X, y=None):
        """Discretizes the input data.

        Parameters
        ----------
        X : array-like, shape [n_samples]
            The array containing continuous features to be discretized.
            Input must be 1d arrays.

        """
        check_is_fitted(self, ["cut_points_"])
        X = column_or_1d(X)
        output = np.searchsorted(self.cut_points_, X)
        return output

