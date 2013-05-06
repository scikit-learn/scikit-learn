import array

import numpy as np
from scipy.sparse import issparse
from scipy.stats import rankdata

from ..base import BaseEstimator, TransformerMixin
from ..utils import atleast2d_or_csc

from .univariate_selection import _BaseFilter


def count_nonzero(X, y=None):
    # In IR/NLP, this is known as Document Frequency (df)
    X = atleast2d_or_csc(X)
    if issparse(X):
        return np.diff(X.indptr)
    return X.astype(bool).sum(axis=1)


def sum_values(X, y=None):
    X = atleast2d_or_csc(X)
    return X.sum(axis=1)


def _calc_percentiles(X, scores):
    ranks = rankdata(scores)  # or use np.unique, bincount, cumsum
    ranks -= ranks.min()
    ranks *= 100. / ranks.max()
    return ranks


class SelectBetween(_BaseFilter):
    def __init__(self, score_func=count_nonzero, minimum=None, maximum=None,
                 scaling=None):
        super(SelectBetween, self).__init__(score_func)
        self.minimum = minimum
        self.maximum = maximum
        self.scaling = scaling

    SCALERS = {
        '%range': lambda X, scores: (scores - scores.min()) * 100 / scores.ptp(),
        '%features': lambda X, scores: scores * 100 / scores.shape[0],
        '%samples': lambda X, scores: scores * 100 / X.shape[0],
        'stds': lambda X, scores: (scores - scores.mean()) / scores.std(),
        'percentile': _calc_percentiles,
        'incrank': lambda X, scores: rankdata(scores),
        'decrank': lambda X, scores: rankdata(-scores),
    }

    def _scale(self, X, scores):
        scaling = self.scaling
        if scaling is None:
            return scores

        if not callable(scaling):
            scaling = self.SCALERS[scaling]
        return scaling(X, scores)

    def fit(self, X, y=None):
        scores = self._scale(self.score_func(X, y))
        self.scaled_scores_ = scores

    def _get_support_mask(self):
        support = np.ones(self.scores_.shape, dtype=bool)
        lo, hi = self.minimum, self.maximum
        if lo is None and hi is None:
            return support

        if lo is not None:
            support &= self.scores_ >= lo
        if hi is not None:
            support &= self.scores_ <= hi
        return support


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
