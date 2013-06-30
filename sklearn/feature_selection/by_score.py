"""Utilities for selecting among features where each is assigned some score

Features can be thus selected by:
    * selecting a limited number of features with the best score
    * applying a lower or upper bound to the score

These bounding parameters may be expressed absolutely or relative to the data
so that they may be meaningful when ported to different data.
"""
# Authors: Joel Nothman, Gilles Louppe
# License: BSD 3 clause

import numbers
from warnings import warn

import numpy as np

from ..externals import six
from ..base import BaseEstimator

from .base import SelectorMixin


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
    if not isinstance(limit, numbers.Integral):
        if limit < 0. or limit > 1.:
            raise ValueError('Float limit must be between 0 and 1')
        limit = int(limit * len(scores))

    n_support = support.sum()
    if n_support <= limit:
        # already meet the limit
        return
    elif limit == 0:
        support[:] = False
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


def mask_by_score(scores, minimum=None, maximum=None, limit=None):
    """Calculate a support mask given a set of feature scores and thresholds

    Parameters
    ----------
    scores : array-like, shape [number of features]
        A real number calculated for each feature representing its importance.

    minimum : float, string or callable, optional
        Where a score is greater than `minimum`, `support` is set to False.
        If a callable, the value is first calculated by applying to `scores`.
        If a string, it is interpreted as the product of '*'-delimited
        expressions, which may be floats or keywords 'mean', 'median', 'min',
        'max', 'sum' or 'size' which evaluate as functions applied to `scores`.

    maximum : float, string or callable, optional
        Where a score is greater than `maximum`, `support` is set to False.
        If a callable, the value is first calculated by applying to `scores`.
        If a string, it is interpreted as the product of '*'-delimited
        expressions, which may be floats or keywords 'mean', 'median', 'min',
        'max', 'sum' or 'size' which evaluate as functions applied to `scores`.

    limit : int, or float, optional
        Limits the number of returned features. If the value is positive, takes
        the `limit` highest-scoring features (after `maximum` and `minimum` are
        applied); if negative, takes the `limit` lowest scoring features.
        A float (from 0 to 1) is treated as a proportion of the number of
        features.

    Returns
    -------
    support : array of bools, shape=[number of features]

    Examples
    --------
    >>> from __future__ import print_function
    >>> print(mask_by_score([5, 3, 2, 4, 1], maximum=4, limit=2))
    ...                                        #doctest: +NORMALIZE_WHITESPACE
    [False True False True False]
    >>> print(mask_by_score([5, 3, 2, 4, 1], minimum='mean'))
    ...                                        #doctest: +NORMALIZE_WHITESPACE
    [ True True False True False]
    >>> print(mask_by_score([5, 3, 2, 4, 1], limit=0.5))
    ...                                        #doctest: +NORMALIZE_WHITESPACE
    [ True False False True False]
    """
    scores = np.asarray(scores)
    support = np.ones(scores.shape, dtype=bool)

    lo = _calc_threshold(scores, minimum)
    hi = _calc_threshold(scores, maximum)
    if lo is not None or hi is not None:
        if lo is not None:
            support &= scores >= lo
        if hi is not None:
            support &= scores <= hi

    if limit is not None:
        _limit_support(support, scores, limit)

    return support


##### Standalone feature selector ######

class SelectByScore(BaseEstimator, SelectorMixin):
    """Select features according to a single score given to each

    Parameters
    ----------
    score_func : callable (X, y) -> array shape=[number of features]
        A function that calculates a real score for each feature representing
        its importance.

    minimum : float, string or callable, optional
        Where a score is greater than `minimum`, `support` is set to False.
        If a callable, the value is first calculated by applying to `scores`.
        If a string, it is interpreted as the product of '*'-delimited
        expressions, which may be floats or keywords 'mean', 'median', 'min',
        'max', 'sum' or 'size' which evaluate as functions applied to `scores`.

    maximum : float, string or callable, optional
        Where a score is greater than `maximum`, `support` is set to False.
        If a callable, the value is first calculated by applying to `scores`.
        If a string, it is interpreted as the product of '*'-delimited
        expressions, which may be floats or keywords 'mean', 'median', 'min',
        'max', 'sum' or 'size' which evaluate as functions applied to `scores`.

    limit : int, or float, optional
        Limits the number of returned features. If the value is positive, takes
        the `limit` highest-scoring features (after `maximum` and `minimum` are
        applied); if negative, takes the `limit` lowest scoring features.
        A float (from 0 to 1) is treated as a proportion of the number of
        features.

    Attributes
    ----------
    scores_ : array, shape=[number of features]
        The scores calculated on the given data.
    """
    def __init__(self, score_func, minimum=None, maximum=None,
                 limit=None):
        self.score_func = score_func
        self.minimum = minimum
        self.maximum = maximum
        self.limit = limit

    def fit(self, X, y=None):
        """Evaluate the score function on samples X with outputs y.

        Parameters
        ----------
        X : array-like, shape=[n_samples, n_features]
            Samples

        y : array-like, shape=[n_features], optional
            Targets
        """
        self.scores_ = self.score_func(X, y)
        return self

    def _get_support_mask(self):
        return mask_by_score(self.scores_,
                             minimum=self.minimum, maximum=self.maximum,
                             limit=self.limit)
