from __future__ import division
from abc import ABCMeta, abstractmethod
import array
from warnings import warn

import numpy as np
from scipy.sparse import issparse
from scipy.stats import rankdata, scoreatpercentile

from ..base import BaseEstimator, TransformerMixin
from ..utils import atleast2d_or_csc, as_float_array
from ..externals import six

from .base import FeatureSelectionMixin


def count_nonzero(X, y=None):
    # In IR/NLP, this is known as Document Frequency (df)
    X = atleast2d_or_csc(X)
    if issparse(X):
        return np.diff(X.indptr)
    return X.astype(bool).sum(axis=1)


def sum_values(X, y=None):
    X = atleast2d_or_csc(X)
    return X.sum(axis=1)


class BaseScaler(object):
    __slots__ = ()
    """
    Scalers, given scores and number of samples, create a closure in which to
    scale some threshold value. If the closure has an attribute
    `scaled_scores`, the threshold will be compared against it.

    These are implemented as classes because closures are not picklable
    """
    def __init__(self, scores, num_samples):
        pass

    def __call__(self, t):
        return t


class _scale_percent_range(BaseScaler):
    def __init__(self, scores, num_samples):
        self.min_ = scores.min()
        self.range_ = scores.ptp()

    def __call__(self, t):
        return t * self.range_ / 100. + self.min_


class _scale_percent_samples(BaseScaler):
    def __init__(self, scores, num_samples):
        self.n_samples = num_samples

    def __call__(self, t):
        return t * self.n_samples / 100.


class _scale_standard_deviations(BaseScaler):
    def __init__(self, scores, num_samples):
        self.mean = np.mean(scores)
        self.std = np.std(scores)

    def __call__(self, t):
        return t * self.std + self.mean


class _scale_percentile(BaseScaler):
    def __init__(self, scores, num_samples):
        self.scores = scores

    def __call__(self, t):
        return scoreatpercentile(self.scores, t)


class _scale_rank(BaseScaler):
    def __init__(self, scores, num_samples):
        self.scaled_scores = rankdata(scores)


class _scale_order(BaseScaler):
    def __init__(self, scores, num_samples):
        self.scaled_scores = np.argsort(scores).argsort()


class SelectByScoreMixin(FeatureSelectionMixin):
    SCALERS = {
        '%range': _scale_percent_range,
        '%samples': _scale_percent_samples,
        'stds': _scale_standard_deviations,
        'percentile': _scale_percentile,
        'incrank': _scale_rank,
        'decrank': lambda scores, n_samples: _scale_rank(-scores, n_samples),
        'incorder': lambda X, scores: scores.argsort().argsort(),
        'decorder': lambda scores, n_samples: _scale_order(-scores, n_samples),
    }

    THRESHOLD_SUBS = {
        'mean': np.mean,
        'median': np.median,
        'min': np.min,
        'max': np.max,
        'sum': np.sum,
        'length': np.size,
    }

    def _make_scaler(self, scaling):
        # TODO: memoize
        if scaling is None:
            scaling = BaseScaler
        elif not callable(scaling):
            scaling = self.SCALERS[scaling]
        return scaling(self.__scores, self.__n_samples)

    def _calc_threshold(self, val):
        if callable(val):
            return val(self.__scores)
        elif isinstance(val, six.string_types):
            res = 1.
            for part in val.split('*'):
                try:
                    res *= float(part)
                except ValueError:
                    try:
                        part = self.THRESHOLD_SUBS[part.strip()]
                    except KeyError:
                        raise ValueError('Unknown reference: %r' % part)
                    else:
                        res *= part(self.__scores)
            return res
        return val

    def _fit(self, scores, X):
        self.__scores = scores
        if X is not None:
            if not issparse(X):
                X = np.asarray(X)
            if len(scores) != X.shape[1]:
                raise ValueError("Scores size differs from number of features")
            self.__n_samples  = X.shape[0]
        else:
            # Hack!
            self.__n_samples = None
        return self

    def _get_support_mask(self, minimum=None, maximum=None, scaling=None, limit=None):
        support = np.ones(self.__scores.shape, dtype=bool)

        lo = self._calc_threshold(minimum)
        hi = self._calc_threshold(maximum)
        if lo is not None or hi is not None:
            scale = self._make_scaler(scaling)
            scaled_scores = getattr(scale, 'scaled_scores', self.__scores)
            if lo is not None:
                support &= scaled_scores >= scale(lo)
            if hi is not None:
                support &= scaled_scores <= scale(hi)

        if limit is not None:
            self._limit_support(support, self.__scores, limit)

        return support

    @staticmethod
    def _limit_support(support, scores, limit):
        """
        Limits the number of True values in support to limit

        Modifies support in-place: if it had more than limit True values,
        it will be left with exactly limit True values.

        A positive limit keeps the highest supported scores, while a negative
        limit keeps the lowest supported scores.

        If limit is a float strictly between 0 and +-1, it will be treated as a
        proportion of the number of scores.
        """
        limit, greatest = abs(limit), limit > 0
        if limit <= 0 or (limit >= 1 and int(limit) != limit):
            raise ValueError(
                'Require 0 < limit < 1 or integer limit >= 1')
        elif limit < 1.:
            limit = int(limit * len(scores))

        n_support = support.sum()
        if n_support <= limit:
            # already meet the limit
            return

        percentile = 100. * min(limit - 1, n_support - 1) / (n_support - 1)
        if greatest:
            percentile = 100. - percentile

        kept_scores = scores[support]
        if issubclass(scores.dtype.type, np.inexact):
            nan_sub = getattr(np.finfo(scores.dtype), 'min' if greatest else 'max')
            kept_scores[np.isnan(kept_scores)] = nan_sub
        cutoff = np.percentile(kept_scores, percentile)

        # take features strictly better than the cutoff
        if greatest:
            support_mask = kept_scores > cutoff
        else:
            support_mask = kept_scores < cutoff

        # Because we chose percentile to exactly match some score, there
        # will always be at least one remaining. Where it matches multiple
        # score, we may need to arbitrarily break the tie.
        n_remaining = limit - support_mask.sum()
        ties = np.where(kept_scores == cutoff)[0]
        if len(ties) > n_remaining:
            warn("Tied features are being arbitrarily split. "
                 "There may be duplicate features, or you used a "
                 "classification score for a regression task.")
        kept_ties = ties[:n_remaining]
        support_mask[kept_ties] = True

        # finally, update support
        support[support] = support_mask


class SelectByScore(BaseEstimator, SelectByScoreMixin):
    def __init__(self, score_func=count_nonzero, minimum=None, maximum=None,
                 scaling=None, limit=None):
        self.score_func = score_func
        self.minimum = minimum
        self.maximum = maximum
        self.scaling = scaling
        self.limit = limit

    def fit(self, X, y=None):
        return super(SelectByScore, self)._fit(self.score_func(X, y), X)

    def _get_support_mask(self):
        return super(SelectByScore, self)._get_support_mask(
            self.minimum, self.maximum, self.scaling, self.limit)


class SelectByName(BaseEstimator, TransformerMixin):
    def __init__(self, transformer, filter, exclude=False):
        self.transformer = transformer
        self.filter = filter
        self.exclude = exclude

    def _get_indices(self):
        filter = self.filter
        if not callable(filter):
            collection = filter
            filter = lambda x: x in collection
        feature_names = self.transformer.get_feature_names()
        indices = array.array("i")
        for i, name in enumerate(feature_names):
            if bool(filter(name)) ^ self.exclude:
                indices.append(i)
        if not indices:
            return np.array([])
        return np.frombuffer(indices, type=np.intc)

    def fit(self, X, y=None):
        self.transformer.fit(X, y)
        self.indices_ = self._get_indices()

    def fit_transform(self, X, y=None):
        X = self.transformer.fit_transform(X, y)
        self.indices_ = self._get_indices()
        return self.transform(X)

    def transform(self, X):
        X = atleast2d_or_csc(X)
        return X[:, self.indices_]

    def get_feature_names(self):
        return np.array(self.transformer.get_feature_names())[self.indices_]


class SelectIndices(BaseEstimator, TransformerMixin):
    def __init__(self, indices):
        self._indices = indices

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        X = atleast2d_or_csc(X)
        return X[:, self.indices_]
