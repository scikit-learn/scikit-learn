# Authors: Luigi De Bianchi <luigi.debianchi@gmail.com>
# License: BSD 3 clause

from __future__ import division

import numpy as np
from scipy import sparse

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.validation import (check_is_fitted, FLOAT_DTYPES)


__all__ = [
    'Smoother',
]

class Smoother(BaseEstimator, TransformerMixin):

    _modes = ["rect"]

    def __init__(self, mode="rect"):
        if mode not in self._modes:
            raise ValueError("Given mode is not supported")
        self._mode = mode

    def _reset(self):
        #TODO: add here all variables of the class
        if hasattr(self, '_mode'):
            del self._mode

    def fit(self, X, y=None):
        self._reset()
        return self.partial_fit(X, y)

    def partial_fit(self, X, y=None):
        if sparse.issparse(X):
            raise TypeError("Smoother does no support sparse input. ")

        X = check_array(X, copy=self.copy, warn_on_dtype=True,
                        estimator=self, dtype=FLOAT_DTYPES)

        # considering a matrix n_samples * n_features
        # the smoothing is applied to every feature

        return self

    def transform(self, X):
        check_is_fitted(self, 'scale_')

        X = check_array(X, copy=self.copy, dtype=FLOAT_DTYPES)



        return X
