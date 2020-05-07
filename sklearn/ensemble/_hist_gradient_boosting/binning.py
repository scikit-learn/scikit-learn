"""
This module contains the BinMapper class.

BinMapper is used for mapping a real-valued dataset into integer-valued bins.
Bin thresholds are computed with the quantiles so that each bin contains
approximately the same number of samples.
"""
# Author: Nicolas Hug

import numpy as np

from ...utils import check_random_state, check_array
from ...base import BaseEstimator, TransformerMixin
from ...utils.validation import check_is_fitted
from ._binning import _map_num_to_bins, _map_cat_to_bins
from .common import X_DTYPE, X_BINNED_DTYPE, ALMOST_INF


def _find_binning_thresholds(data, max_bins, is_categorical=None):
    """Extract feature-wise quantiles from numerical data.

    Missing values are ignored for finding the thresholds.

    Parameters
    ----------
    data : array-like, shape (n_samples, n_features)
        The data to bin.
    max_bins: int
        The maximum number of bins to use for non-missing values. If for a
        given feature the number of unique values is less than ``max_bins``,
        then those unique values will be used to compute the bin thresholds,
        instead of the quantiles
    is_categorical : ndarray of bool or None
        Indicates categorical features

    Return
    ------
    binning_thresholds: list of arrays
        For each feature, stores the increasing numeric values that can
        be used to separate the bins. Thus ``len(binning_thresholds) ==
        n_features``.
    """
    binning_thresholds = []
    for f_idx in range(data.shape[1]):
        # categorical feature
        if is_categorical is not None and is_categorical[f_idx]:
            binning_thresholds.append(None)
            continue

        col_data = data[:, f_idx]
        # ignore missing values when computing bin thresholds
        missing_mask = np.isnan(col_data)
        if missing_mask.any():
            col_data = col_data[~missing_mask]
        col_data = np.ascontiguousarray(col_data, dtype=X_DTYPE)
        distinct_values = np.unique(col_data)
        if len(distinct_values) <= max_bins:
            midpoints = distinct_values[:-1] + distinct_values[1:]
            midpoints *= .5
        else:
            # We sort again the data in this case. We could compute
            # approximate midpoint percentiles using the output of
            # np.unique(col_data, return_counts) instead but this is more
            # work and the performance benefit will be limited because we
            # work on a fixed-size subsample of the full data.
            percentiles = np.linspace(0, 100, num=max_bins + 1)
            percentiles = percentiles[1:-1]
            midpoints = np.percentile(col_data, percentiles,
                                      interpolation='midpoint').astype(X_DTYPE)
            assert midpoints.shape[0] == max_bins - 1

        # We avoid having +inf thresholds: +inf thresholds are only allowed in
        # a "split on nan" situation.
        np.clip(midpoints, a_min=None, a_max=ALMOST_INF, out=midpoints)

        binning_thresholds.append(midpoints)

    return binning_thresholds


def _find_categories(data, max_bins, is_categorical):
    """Extract feature-wise categories from categorical data

    Missing values and negative values are ignored. They will be considered
    missing when ``_encode_categories`` is called.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        The data to bin.
    max_bins: int
        The maximum number of bins to use for non-missing values. If for a
        given feature the number of unique values is less than ``max_bins``,
        then those unique values will be used to compute the bin thresholds,
        instead of the quantiles
    is_categorical : ndarray of bool of shape (n_features,)
        Indicates categorical features.

    Return
    ------
    bin_categories : list of arrays or None
        For each categorical feature, this gives the categories corresponding
        to each bin.
    """
    data = data[:, is_categorical]
    bin_categories = []
    for f_idx in range(data.shape[1]):
        col_data = data[:, f_idx]

        categories, counts = np.unique(col_data, return_counts=True)

        # sort by highest count
        sorted_idx = np.argsort(-counts)
        categories = categories[sorted_idx]

        # nans and negative values will be considered missing
        missing = np.isnan(categories)
        negative = categories < 0
        both = missing | negative
        if both.any():
            categories = categories[~both]

        # keep at most max_bins categories
        # needs to be sorted, because `_encode_categories` will use
        # np.searchsorted for encoding
        bin_categories.append(np.sort(categories[:max_bins]))

    return bin_categories


class _BinMapper(TransformerMixin, BaseEstimator):
    """Transformer that maps a dataset into integer-valued bins.

    The bins are created in a feature-wise fashion, using quantiles so that
    each bins contains approximately the same number of samples.

    For large datasets, quantiles are computed on a subset of the data to
    speed-up the binning, but the quantiles should remain stable.

    Features with a small number of values may be binned into less than
    ``n_bins`` bins. The last bin (at index ``n_bins - 1``) is always reserved
    for missing values.

    Parameters
    ----------
    n_bins : int, optional (default=256)
        The maximum number of bins to use (including the bin for missing
        values). Non-missing values are binned on ``max_bins = n_bins - 1``
        bins. The last bin is always reserved for missing values. If for a
        given feature the number of unique values is less than ``max_bins``,
        then those unique values will be used to compute the bin thresholds,
        instead of the quantiles. For categorical features indicated by
        ``is_categorical``, the docstring for ``is_categorical`` details on
        this procedure.
    subsample : int or None, optional (default=2e5)
        If ``n_samples > subsample``, then ``sub_samples`` samples will be
        randomly chosen to compute the quantiles. If ``None``, the whole data
        is used.
    is_categorical : ndarray of bool of shape (n_features,), default=None
        Indicates categorical features. If the cardinality of a categorical
        feature is greater than ``n_bins``, then the ``n_bins`` most frequent
        categories are kept. The infrequent categories will be consider
        missing. During ``transform`` time, unknown categories will also be
        considered missing.
    random_state: int, RandomState instance or None
        Pseudo-random number generator to control the random sub-sampling.
        Pass an int for reproducible output across multiple
        function calls.
        See :term: `Glossary <random_state>`.

    Attributes
    ----------
    bin_thresholds_ : list of arrays
        For each feature, gives the real-valued bin threhsolds. There are
        ``max_bins - 1`` thresholds, where ``max_bins = n_bins - 1`` is the
        number of bins used for non-missing values.
        If ``bin_thresholds_`` corresponds to a categorical feature, then
        ``bin_thresholds_[categorical_index]`` is None.
    bin_categories_ : list of arrays or empty list
        For each categorical feature, this gives the categories corresponding
        to each bin.
    n_bins_non_missing_ : array of uint32
        For each feature, gives the number of bins actually used for
        non-missing values. For features with a lot of unique values, this is
        equal to ``n_bins - 1``.
    missing_values_bin_idx_ : uint8
        The index of the bin where missing values are mapped. This is a
        constant across all features. This corresponds to the last bin, and
        it is always equal to ``n_bins - 1``. Note that if ``n_bins_missing_``
        is less than ``n_bins - 1`` for a given feature, then there are
        empty (and unused) bins.
    """
    def __init__(self, n_bins=256, subsample=int(2e5), is_categorical=None,
                 random_state=None):
        self.n_bins = n_bins
        self.subsample = subsample
        self.is_categorical = is_categorical
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit data X by computing the binning thresholds.

        The last bin is reserved for missing values, whether missing values
        are present in the data or not.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The data to bin.
        y: None
            Ignored.

        Returns
        -------
        self : object
        """
        if not (3 <= self.n_bins <= 256):
            # min is 3: at least 2 distinct bins and a missing values bin
            raise ValueError('n_bins={} should be no smaller than 3 '
                             'and no larger than 256.'.format(self.n_bins))

        X = check_array(X, dtype=[X_DTYPE], force_all_finite=False)
        max_bins = self.n_bins - 1

        rng = check_random_state(self.random_state)
        if self.subsample is not None and X.shape[0] > self.subsample:
            subset = rng.choice(X.shape[0], self.subsample, replace=False)
            X = X.take(subset, axis=0)

        self.bin_thresholds_ = _find_binning_thresholds(
            X, max_bins, is_categorical=self.is_categorical)

        if self.is_categorical is not None and np.any(self.is_categorical):
            self.bin_categories_ = _find_categories(
                X, max_bins, is_categorical=self.is_categorical)
        else:
            self.bin_categories_ = []

        n_bins_non_missing = []

        if self.bin_categories_:
            categorical_indices = np.flatnonzero(self.is_categorical)
            cat_idx_to_bin = dict(zip(categorical_indices,
                                      self.bin_categories_))

        for i, thresholds in enumerate(self.bin_thresholds_):
            if self.bin_categories_ and i in cat_idx_to_bin:
                # category
                n_bins_non_missing.append(cat_idx_to_bin[i].shape[0])
            else:
                # numerical
                n_bins_non_missing.append(thresholds.shape[0] + 1)

        self.n_bins_non_missing_ = np.array(n_bins_non_missing,
                                            dtype=np.uint32)

        self.missing_values_bin_idx_ = self.n_bins - 1

        return self

    def transform(self, X):
        """Bin data X.

        Missing values will be mapped to the last bin.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The data to bin.

        Returns
        -------
        X_binned : array-like, shape (n_samples, n_features)
            The binned data (fortran-aligned).
        """
        X = check_array(X, dtype=[X_DTYPE], force_all_finite=False)
        check_is_fitted(self)
        if X.shape[1] != self.n_bins_non_missing_.shape[0]:
            raise ValueError(
                'This estimator was fitted with {} features but {} got passed '
                'to transform()'.format(self.n_bins_non_missing_.shape[0],
                                        X.shape[1])
            )

        binned = np.zeros_like(X, dtype=X_BINNED_DTYPE, order='F')
        _map_num_to_bins(X, self.bin_thresholds_, self.missing_values_bin_idx_,
                         binned)
        if self.bin_categories_:
            categorical_indices = np.flatnonzero(self.is_categorical)
            _map_cat_to_bins(X, categorical_indices,
                             self.bin_categories_,
                             self.missing_values_bin_idx_, binned)
        return binned

    def transform_categories_only(self, X):
        """Bin data in X that is categorical.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The data to bin.

        Returns
        -------
        X_binned : array-like, shape (n_samples, n_categorical_features)
            The binned data (F-aligned)
        """
        check_is_fitted(self)
        if not self.bin_categories_:
            raise ValueError("transform_categories_only can only be set when "
                             "there are categorical features in fit")
        X = check_array(X[:, self.is_categorical], dtype=[X_DTYPE],
                        force_all_finite=False)

        n_samples = X.shape[0]
        n_categories = len(self.bin_categories_)
        binned = np.zeros((n_samples, n_categories), dtype=X_BINNED_DTYPE,
                          order='F')
        _map_cat_to_bins(X, np.arange(n_categories),
                         self.bin_categories_,
                         self.missing_values_bin_idx_, binned)
        return binned
