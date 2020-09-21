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
from ._binning import _map_to_bins
from .common import X_DTYPE, X_BINNED_DTYPE, ALMOST_INF, X_BITSET_INNER_DTYPE
from ._bitset import set_bitset_memoryview


def _find_binning_threshold(col_data, max_bins):
    """Extract feature-wise quantiles from numerical data.

    Missing values are ignored for finding the thresholds.

    Parameters
    ----------
    col_data : array-like, shape (n_features,)
        The numerical feature to bin.
    max_bins: int
        The maximum number of bins to use for non-missing values. If for a
        given feature the number of unique values is less than ``max_bins``,
        then those unique values will be used to compute the bin thresholds,
        instead of the quantiles

    Return
    ------
    binning_thresholds: ndarray
        For each feature, stores the increasing numeric values that can
        be used to separate the bins. Thus ``len(binning_thresholds) ==
        n_features``.
    """
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
    return midpoints


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
    n_bins : int, default=256
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
        categories are kept. During ``transform`` time, unknown categories will
        be considered missing.
    random_state: int, RandomState instance or None, default=None
        Pseudo-random number generator to control the random sub-sampling.
        Pass an int for reproducible output across multiple
        function calls.
        See :term: `Glossary <random_state>`.

    Attributes
    ----------
    bin_thresholds_ : list of ndarray
        For each feature, each array indicates how to map a feature into a
        binned feature. The semantic and size depends on the nature of the
        feature:
        - for real-valued features, the array corresponds to the real-valued
          bin thresholds (the upper bound of each bin). There are ``max_bins
          - 1`` thresholds, where ``max_bins = n_bins - 1`` is the number of
          bins used for non-missing values.
        - for categorical features, the array is a map from a binned category
          value to the raw category value. The size of the array is equal to
          ``min(max_bins, category_cardinality)`` where we ignore missing
          values in the cardinality.
    n_bins_non_missing_ : ndarray, dtype=np.uint32
        For each feature, gives the number of bins actually used for
        non-missing values. For features with a lot of unique values, this is
        equal to ``n_bins - 1``.
    is_categorical_ : ndarray of shape (n_features,), dtype=np.uint8
        Indicator for categorical features.
    missing_values_bin_idx_ : np.uint8
        The index of the bin where missing values are mapped. This is a
        constant across all features. This corresponds to the last bin, and
        it is always equal to ``n_bins - 1``. Note that if ``n_bins_missing_``
        is less than ``n_bins - 1`` for a given feature, then there are
        empty (and unused) bins.
    """
    def __init__(self, n_bins=256, subsample=int(2e5), is_categorical=None,
                 known_categories=None, random_state=None):
        self.n_bins = n_bins
        self.subsample = subsample
        self.is_categorical = is_categorical
        self.known_categories = known_categories
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit data X by computing the binning thresholds.

        The last bin is reserved for missing values, whether missing values
        are present in the data or not.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
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

        if self.is_categorical is None:
            self.is_categorical_ = np.zeros(X.shape[1], dtype=np.uint8)
        else:
            self.is_categorical_ = np.asarray(self.is_categorical,
                                              dtype=np.uint8)

        known_categories = self.known_categories
        if known_categories is None:
            known_categories = [None] * X.shape[1]

        for f_idx, (is_categorical, cats) in enumerate(zip(
                self.is_categorical_,
                known_categories)):
            if is_categorical and cats is None:
                raise ValueError(
                    f"Known categories for feature {f_idx} must be provided."
                )
            if cats is not None and not is_categorical:
                raise ValueError(
                    f"Feature {f_idx} isn't marked as a categorical feature, "
                    f"but categories where passed."
                )

        self.missing_values_bin_idx_ = self.n_bins - 1

        bin_thresholds = []
        n_bins_non_missing = []

        for f_idx in range(X.shape[1]):
            if not self.is_categorical_[f_idx]:
                bins = _find_binning_threshold(X[:, f_idx], max_bins)
                n_bins_non_missing.append(bins.shape[0] + 1)
            else:
                bins = known_categories[f_idx]
                n_bins_non_missing.append(bins.shape[0])

            bin_thresholds.append(bins)

        self.bin_thresholds_ = bin_thresholds
        self.n_bins_non_missing_ = np.array(n_bins_non_missing,
                                            dtype=np.uint32)
        return self

    def transform(self, X):
        """Bin data X.

        Missing values will be mapped to the last bin.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to bin.

        Returns
        -------
        X_binned : array-like of shape (n_samples, n_features)
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
        _map_to_bins(X, self.bin_thresholds_, self.missing_values_bin_idx_,
                     self.is_categorical_, binned)
        return binned

    def make_known_categories(self):
        """Create a mapping from original feature indices to known categorical
        indices and a bitset for known categories.

        Returns
        -------
        result: dict
            Dictionary with items:
            - `known_cat_bitset` : ndarray of shape (n_categorical, 8)
                Bitset for known categories.
            - `orig_feat_to_known_cat_idx`: ndarray of shape (n_features,)
                Maps from original feature indies to known categorical
                indices.
        """

        n_features = len(self.bin_thresholds_)
        categorical_indices = np.flatnonzero(self.is_categorical_)
        orig_feat_to_known_cats_idx = np.zeros(n_features, dtype=np.uint8)

        orig_feat_to_known_cats_idx[categorical_indices] = np.arange(
            categorical_indices.size, dtype=np.uint8)

        known_cat_bitset = np.zeros((categorical_indices.size, 8),
                                    dtype=X_BITSET_INNER_DTYPE)
        for idx, cat_idx in enumerate(categorical_indices):
            for num in self.bin_thresholds_[cat_idx]:
                set_bitset_memoryview(known_cat_bitset[idx], num)

        return {
            'known_cat_bitset': known_cat_bitset,
            'orig_feat_to_known_cats_idx': orig_feat_to_known_cats_idx
        }
