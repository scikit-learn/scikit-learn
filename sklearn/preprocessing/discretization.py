# Author: Henry Lin <hlin117@gmail.com>

from __future__ import division

from itertools import izip
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_array, check_is_fitted

class KBinsDiscretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    n_bins : int (default=2)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.
        Raises ValueError if n_bins < 2.

    Attributes
    ----------
    min_ : float
        The minimum value of the input data.

    max_ : float
        The maximum value of the input data.

    n_features_ : int
        The number of features from the original dataset.

    cut_points_ : array of shape (n_bins - 1, n_features_)
        Each column of cut_points_ is a set of cut points which describe
        how the values of the corresponding feature are binned. For a pair
        of adjacent values in a column, a, b, values between [a, b) will
        have the same discretization value.

    Example
    -------
    >>> X = [[-2, 1, -4,   -1], \
             [-1, 2, -3, -0.5], \
             [ 0, 3, -2,  0.5], \
             [ 1, 4, -1,    2]]
    >>> dis = KBinsDiscretizer(n_bins=3)
    >>> dis.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)
    >>> dis.cut_points_
    array([[-1.,  2., -3.,  0.],
           [ 0.,  3., -2.,  1.]])
    >>> dis.transform(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[-2., 1., -4., -1.],
           [-1., 2., -3., -1.],
           [ 0., 3., -2.,  0.],
           [ 0., 3., -2.,  1.]])
    """

    sparse_formats = ['csr', 'csc']


    def __init__(self, n_bins=2):
        self.copy = copy
        self.n_bins = n_bins

        # Attributes
        self.min_ = None
        self.max_ = None
        self.n_features_ = None


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
                             "of bins. Received {0}, expected at least 2."
                             .format(self.n_bins))

        X = check_array(X, accept_sparse=self.sparse_formats, dtype=float,
                        ensure_2d=False)
        self.n_features_ = X.shape[1]

        self.min_ = X.min(axis=0)
        self.max_ = X.max(axis=0)
        return self


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
        check_is_fitted(self, ["min_", "max_"])
        X = check_array(X, accept_sparse=self.sparse_formats, dtype=float,
                        ensure_2d=False)

        if X.shape[1] != self.n_features_:
            raise ValueError("Transformed array does not have currect number "
                             "of features. Expecting {0}, received {1}."
                             .format(self.n_features_, X.shape[1]))

        # TODO: Make copy parameter?
        X = X.copy()

        # Clip outliers from data (TODO)
        #X[X > self.max_] = self.max_
        #X[X < self.min_] = self.min_

        return self._transform(X)


    def _transform(self, X):
        """Discretizes the continuous features into k integers.

        The process used below preserves the sparsity of the data by
        transforming the data, and rounding the data down to the
        nearest integer. These integers represent the bins sought.

        The process is as follows.
        1. Let X_trans = X * n_bins / (X.max() - X.min()).
            Observe that X_trans is a dataset of width k.
        2. epsilon = abs(X_trans.min() - ceil(X_trans.min())).
            episilon < 1 is the distance of X.min() to ceil(X.min())
        3. X_trans[X_trans != 0] += epsilon
            Shifts the data. X_trans.min() is now an integer. Note this
            operation can be done on sparse matrices.
        4. output = floor(X_trans)
            The output now contains integers which represent the
            discretized values. Observe that even if we added epsilon to
            zero values in step 3, these zeros will still be discretized
            to zero.
        """
        # Rescale data
        X_trans = X / (self.max_ - self.min_)
        X_trans *= self.n_bins

        # Shift data
        epsilon = np.abs(X_trans.min() - np.ceil(X_trans.min()))
        X_trans[X_trans != 0] += epsilon

        # Generate output
        discretized = np.floor(X_trans)

        # Addressing corner case issues
        for col in discretized.T:
            values = np.unique(col)
            col_min = values[0]
            col_max = values[-1]
            if col_max - col_min > self.n_bins - 1:
                col_max = values[-1]
                col[col == col_max] = col_max - 1

        return discretized


    @property
    def cut_points_(self):
        check_is_fitted(self, ["min_", "max_"])

        n_bins = self.n_bins
        mins_ = self.min_
        maxes_ = self.max_

        spacings = (maxes_ - mins_) / n_bins
        cut_points = [np.linspace(min_ + spacing, max_ - spacing, num=n_bins-1)
                      for spacing, min_, max_ in izip(spacings, mins_, maxes_)]
        return np.array(cut_points).T
