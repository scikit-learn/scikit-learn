# Title          : discretization.py
# Author         : Henry Lin <hlin117@gmail.com>
# License        : BSD 3 clause
#==============================================================================

from __future__ import division
import numpy as np
from six.moves import xrange
from ..base import BaseEstimator, TransformerMixin
from ..utils import column_or_1d
from sklearn.utils.validation import check_is_fitted

class FixedWidthDiscretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    numBins : int (default=2)
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

    Example
    -------
    >>> X = [1, 2, 3, 4, 5, 6]
    >>> discretizer = FixedWidthDiscretizer(nbins=3)
    >>> discretizer.fit(X) # DOCTEST +DONT_ACCEPT_BLANKLINE
    >>> discretizer.cut_points_
    array([2., 4.]...)
    >>> discretizer.transform(X)
    array([0, 0, 1, 1, 2, 2])
    """

    def __init__(self, numBins=2):
        if numBins < 2:
            raise ValueError("FixedWidthDiscretizer received an invalid number "
                             "of bins")
        self.numBins = numBins

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
        self.spacing_ = (self.max_ - self.min_) / self.numBins
        cut_points = (self.min_ + self.spacing_ * i for i in xrange(1, self.numBins))
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

