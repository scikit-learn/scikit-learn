"""Utilities for selecting among features where each is assigned some score

Features can be thus selected by:
    * selecting a limited number of features with the best score
    * applying a lower or upper bound to the score

These bounding parameters may be expressed absolutely or relative to the data
so that they may be meaningful when ported to different data.  This motivates
both flexible specification of the bounds, and a means to reinterpret them by
transforming them into a different space using "scalers".
"""

from warnings import warn

import numpy as np
from scipy.sparse import issparse
from scipy.stats import rankdata, scoreatpercentile

from ..externals import six
from ..utils import atleast2d_or_csc
from ..base import BaseEstimator

from .base import FeatureSelectionMixin


###### Simple scores ######

def count_nonzero(X, y=None):
    # In IR/NLP, this is known as Document Frequency (df)
    X = atleast2d_or_csc(X)
    if issparse(X):
        return np.diff(X.indptr)
    return X.astype(bool).sum(axis=1)


def sum_values(X, y=None):
    X = atleast2d_or_csc(X)
    return X.sum(axis=1)


###### Scalers ######

class BaseScaler(object):
    """An extensible implementation of null scaling

    Scalers, given scores and number of samples, create a closure in which to
    scale some threshold value. If the closure has an attribute
    `scaled_scores`, the threshold will be compared against it.

    These are implemented as classes because closures are not picklable
    """

    __slots__ = ()

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


###### Strings interpreted by mask_by_score ######

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


###### mask_by_score and its helpers ######

def _calc_threshold(scores, val):
    if callable(val):
        return val(scores)
    elif isinstance(val, six.string_types):
        res = 1.
        for part in val.split('*'):
            try:
                res *= float(part)
            except ValueError:
                try:
                    part = THRESHOLD_SUBS[part.strip()]
                except KeyError:
                    raise ValueError('Unknown reference: %r' % part)
                else:
                    res *= part(scores)
        return res
    return val


def _make_scaler(scores, n_samples, scaling):
    # TODO: memoize
    if scaling is None:
        scaling = BaseScaler
    elif not callable(scaling):
        scaling = SCALERS[scaling]
    return scaling(scores, n_samples)


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


def mask_by_score(scores, X_shape=None, minimum=None, maximum=None,
                  scaling=None, limit=None):
    """Calculate a support mask given a set of feature scores and thresholds

    Parameters
    ----------
    scores : array-like, shape [number of features]
        A real number calculated for each feature representing its importance.

    X_shape : pair of integers, optional
        This is the shape of the training data, being the number of samples
        and the number of features.  If present the number of samples may be
        used in scaling, and the number of features is used for validation.

    minimum : float, string or callable, optional
        If specified, the value is scaled accoring to `scaling`, and where a
        score is greater than the scaled minimum, `support` is set to False.

        If a callable, the value is first calculated by applying to `scores`.
        If a string, it is interpreted as the product of '*'-delimited
        expressions, which may be floats or keywords 'mean', 'median', 'min',
        'max', 'sum' or 'size' which evaluate as functions applied to `scores`.

    maximum : float, string or callable, optional
        If specified, the value is scaled accoring to `scaling`, and where a
        score is greater than the scaled maximum, `support` is set to False.

        If a callable, the value is first calculated by applying to `scores`.
        If a string, it is interpreted as the product of '*'-delimited
        expressions, which may be floats or keywords 'mean', 'median', 'min',
        'max', 'sum' or 'size' which evaluate as functions applied to `scores`.

    scaling : callable or string, optional
        As a callable, `minimum` is replaced with `scaling(scores, X_shape[1])(minimum)`, and
        similar for `maximum`.  If `scaling` has an attribute `scaled_scores`,
        `minimum` and `maximum` are compared to that value rather than
        `scores`.  A string value will apply a predefined scaler:
            - 

    Returns
    -------
    support : 
    """
    scores = np.asarray(scores)
    if X_shape is not None:
        if len(scores) != X_shape[1]:
            raise ValueError("Scores size differs from number of features")
        n_samples = X_shape[0]
    else:
        n_samples = None
    support = np.ones(scores.shape, dtype=bool)

    lo = _calc_threshold(scores, minimum)
    hi = _calc_threshold(scores, maximum)
    if lo is not None or hi is not None:
        scale = _make_scaler(scores, n_samples, scaling)
        scaled_scores = getattr(scale, 'scaled_scores', scores)
        if lo is not None:
            support &= scaled_scores >= scale(lo)
        if hi is not None:
            support &= scaled_scores <= scale(hi)

    if limit is not None:
        _limit_support(support, scores, limit)

    return support


##### Standalone feature selector ######

class SelectByScore(BaseEstimator, FeatureSelectionMixin):
    def __init__(self, score_func=count_nonzero, minimum=None, maximum=None,
                 scaling=None, limit=None):
        self.score_func = score_func
        self.minimum = minimum
        self.maximum = maximum
        self.scaling = scaling
        self.limit = limit

    def fit(self, X, y=None):
        self.scores_ = self.score_func(X, y)
        if not issparse(X):
            X = np.asarray(X)
        self._X_shape = X.shape
        return self

    def _get_support_mask(self):
        return mask_by_score(self.scores_, self._X_shape,
                             minimum=self.minimum, maximum=self.maximum,
                             scaling=self.scaling, limit=self.limit)

