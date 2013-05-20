from __future__ import division
import array

import numpy as np

from ..utils import atleast2d_or_csc
from ..base import BaseEstimator, TransformerMixin


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
