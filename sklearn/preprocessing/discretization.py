# Author: Henry Lin <hlin117@gmail.com>

from __future__ import division
import numpy as np
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.sparsefuncs import csr_csc_min_axis0, csr_csc_max_axis0
from _discretization import binsearch
from ..externals.six.moves import xrange

__all__ = [
    "Discretizer"
]


class Discretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    n_bins : int (default=2)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.
        Raises ValueError if n_bins < 2.

    categorical_features : int, array-like, or None (default=None)
        Specifies the indices of categorical columns which are not to
        be discretized. If None, assumes that all columns are continuous
        features.

    Attributes
    ----------
    min_ : float
        The minimum value of the input data.

    max_ : float
        The maximum value of the input data.

    cut_points_ : array, shape [n_bins - 1, n_continuous_features]
        Contains the boundaries for which the data lies. Each interval
        has a closed left boundary, and an open right boundary.

        Given a feature, the width of each interval is given by
        (max - min) / n_bins.

    zero_intervals_ : list of tuples of length n_continuous_features
        A list of 2-tuples that represents the intervals for which a number
        would be discretized to zero. If an interval is given by (a, b),
        then a value between [a, b) will be discretized to zero.

    searched_points_ : array, shape [numBins - 2, n_continuous_features]
        An array of cut points used for discretization. This array is empty
        if n_bins = 2.

    n_features_ : int
        The number of features from the original dataset.

    n_continuous_features_ : int
        The number of continuous features.

    continuous_features_ : list
        Contains the indices of continuous columns in the dataset.
        This list is sorted.

    Example
    -------
    X has three examples, with X[:, 2] as categorical

    >>> X = [[-3, 1, 0, 5  ], \
             [-2, 7, 8, 4.5], \
             [3,  3, 1, 4  ]]
    >>> from sklearn.preprocessing import Discretizer
    >>> discretizer = Discretizer(n_bins=4, categorical_features=[2])
    >>> discretizer.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Discretizer(...)
    >>> discretizer.cut_points_ # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[-1.5 ,  2.5 ,  4.25],
           [ 0.  ,  4.  ,  4.5 ],
           [ 1.5 ,  5.5 ,  4.75]])

    >>> discretizer.zero_intervals_
    [(0.0, 1.5), (-inf, 2.5), (-inf, 4.25)]

    >>> discretizer.searched_points_ # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[-1.5 ,  4.  ,  4.5 ],
           [ 1.5 ,  5.5 ,  4.75]])

    Transforming X will move the categorical features to the last indices

    >>> discretizer.transform(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[ 1.,  0.,  3.,  0.],
           [ 1.,  3.,  2.,  8.],
           [ 3.,  1.,  0.,  1.]])
    """
    sparse_formats = ['csc'] #['csr', 'csc']

    def __init__(self, n_bins=2, categorical_features=None):
        self.n_bins = n_bins
        self.categorical_features = categorical_features

        # Attributes
        self.min_ = None
        self.max_ = None
        self.n_features_ = None
        self.continuous_features_ = categorical_features
        self.searched_points_ = None
        self.zero_intervals_ = None

    def _set_continuous_features(self):
        """Sets a boolean array that determines which columns are
        continuous features.
        """

        # Case when all features are continuous
        if self.categorical_features is None:
            self.continuous_features_ = range(self.n_features_)
            return

        if type(self.categorical_features) is int:
            categorical_features = [self.categorical_features]
        else:
            categorical_features = self.categorical_features

        if len(categorical_features) > self.n_features_:
            raise ValueError("The number of categorical indices is more than "
                             "the number of features in the dataset.")

        self.continuous_features_ = [i for i in xrange(self.n_features_)
                                     if i not in set(categorical_features)]

        # Checks if there are duplicate columns from self.categorical_features,
        # such as [-1, 0, 0, self.n_features - 1]
        if self.n_features_ - len(categorical_features) \
                != len(self.continuous_features_):
            raise ValueError("Duplicate indices detected from input "
                             "categorical indices. Input was: {0}" \
                             .format(self.continuous_features_))

    def fit(self, X, y=None):
        """Finds the intervals of interest from the input data.

        Parameters
        ----------
        X : {array-like, sparse-matrix}, shape [n_samples, n_features]
            The input data containing continuous features to be
            discretized. Categorical columns are indicated by
            categorical_columns from the constructor.

        Returns
        -------
        self
        """

        if self.n_bins < 2:
            raise ValueError("Discretizer received an invalid number "
                             "of bins. Received {0}, expected at least 2"
                             .format(self.n_bins))

        X = check_array(X, accept_sparse=self.sparse_formats)
        self.n_features_ = X.shape[1]

        # Set the indices of continuous features
        self._set_continuous_features()

        continuous = X[:, self.continuous_features_]

        if sp.issparse(X):
            self.min_ = csr_csc_min_axis0(continuous)
            self.max_ = csr_csc_max_axis0(continuous)
        else:
            self.min_ = continuous.min(axis=0)
            self.max_ = continuous.max(axis=0)

        searched_points = list()
        zero_intervals = list()

        for min_, max_ in zip(self.min_, self.max_):

            points = np.linspace(min_, max_, num=self.n_bins, endpoint=False)[1:]

            # Get index of where zero goes. Omit this index in
            # the rebuilt array
            # TODO: Watch out for when there is only two intervals
            zero_index = np.searchsorted(points, 0, side='right')

            # Case when all values are positive
            if zero_index == 0:
                zero_int = (-np.inf, points[zero_index])
                searched = points[1:]

            # Case when all values are negative
            elif zero_index == len(points):
                zero_int = (points[zero_index - 1], np.inf)
                searched = points[:-1]

            else:
                zero_int = (points[zero_index - 1], points[zero_index])
                searched = np.hstack((points[:zero_index - 1], points[zero_index:]))
            zero_intervals.append(zero_int)

            searched_points.append(searched.reshape(-1, 1))

        self.searched_points_ = np.hstack(searched_points)
        self.zero_intervals_ = zero_intervals
        return self

    @property
    def n_continuous_features_(self):
        if not self.continuous_features_:
            return None
        return len(self.continuous_features_)

    @property
    def cut_points_(self):
        if not self.zero_intervals_:
            return None
        zero_intervals = self.zero_intervals_
        searched_points = self.searched_points_

        cut_points = list()
        for (lower, upper), col in zip(zero_intervals, searched_points.T):

            if lower == -np.inf:
                cut_column = np.insert(col, 0, upper)
            elif upper == np.inf:
                cut_column = np.append(col, lower)
            else:
                zero_index = np.searchsorted(col, lower)
                cut_column = np.insert(col, zero_index, lower)
            cut_points.append(cut_column.reshape(-1, 1))
        cut_points = np.hstack(cut_points)
        return cut_points

    @property
    def cat_features_(self):
        """Returns a list of categorical features, in sorted order.
        """
        if type(self.categorical_features) is int:
            return [self.categorical_features]
        else:
            return sorted(self.categorical_features)


    def _transform_sparse(self, X):
        """Helper function to transform sp.csc_matrices.
        """

        def get_csc_columns(X):
            n_features = X.shape[1]
            col_ptr = X.indptr
            for i in xrange(n_features):
                start = col_ptr[i]
                end = col_ptr[i + 1]
                yield start, end, X.data[start:end]

        searched = self.searched_points_
        z_intervals = self.zero_intervals_

        continuous = X[:, self.continuous_features_]

        # Indices get mixed up using the selector above; the data needs
        # to be sorted for sparse discretization.
        continuous.sort_indices()
        output = continuous.copy()

        if searched.size == 0:
            col_groups = zip(get_csc_columns(continuous), z_intervals)
            for (start, end, cont), (lower, upper) in col_groups:
                zero_mask = np.logical_and(lower <= cont, cont < upper)
                view = output.data[start:end]
                view[zero_mask] = 0
                view[~zero_mask] = 1
        else:
            col_groups = zip(get_csc_columns(continuous), z_intervals, searched.T)
            for (start, end, cont), z_int, points in col_groups:
                dis = binsearch(cont, z_int, points)
                output.data[start:end] = dis

        if self.n_continuous_features_ == self.n_features_:
            return output

        cat_features = self.cat_features_
        categorical = X[:, cat_features]
        return sp.hstack((output, categorical))

    def transform(self, X, y=None):
        """Discretizes the input data.

        Parameters
        ----------
        X : {array-like, sparse-matrix}, shape [n_samples, n_features]
            The input data containing continuous features to be
            discretized. Categorical columns are indicated by
            categorical_columns from the constructor.

        Returns
        -------
        T : array-like, shape [n_samples, n_features]
            Array containing discretized features alongside the
            categorical features. Categorical features are represented
            as the last features of the output.

        """
        check_is_fitted(self, ["cut_points_", "continuous_features_"])
        X = check_array(X, accept_sparse=self.sparse_formats)
        if X.shape[1] != self.n_features_:
            raise ValueError("Transformed array does not have currect number "
                             "of features. Expecting {0}, received {1}"
                             .format(self.n_features_, X.shape[1]))

        if isinstance(X, sp.csc_matrix):
            return self._transform_sparse(X)

        continuous = X[:, self.continuous_features_]
        discretized = list()
        z_intervals = self.zero_intervals_
        searched = self.searched_points_

        # Case when n_bins == 2
        if searched.size == 0:
            for cont, (lower, upper) in zip(continuous.T, z_intervals):
                dis = np.zeros(len(cont))
                zero_mask = np.logical_and(lower <= cont, cont < upper)
                dis[~zero_mask] = 1
                discretized.append(dis.reshape(-1, 1))
        else:
            for cont, z_int, points in zip(continuous.T, z_intervals, searched.T):
                dis = binsearch(cont, z_int, points)
                discretized.append(dis.reshape(-1, 1))

        discretized = np.hstack(discretized)
        if self.n_continuous_features_ == self.n_features_:
            return discretized

        cat_features = self.cat_features_
        categorical = X[:, cat_features]
        return np.hstack((discretized, categorical))

