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

    def __init__(self, mode="rect", size=3):
        if mode not in self._modes:
            raise ValueError("Given mode is not supported.")

        if size%2 == 0:
            raise ValueError("Size must be odd.")

        self._mode = mode
        self._size = size

    def _reset(self):
        if hasattr(self, '_mode'):
            del self._mode
            del self._size

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

        wiskersLen = int((self._size - 1) / 2)
        sampleSize = X.shape[0]

        if self._mode == "rect":
            for row in range(wiskersLen, sampleSize-wiskersLen):
                #print("row: " + str(row) + " \t ")
                #print(X[row, :])
                average = self._sizeSum(X, row - wiskersLen, row + wiskersLen + 1) / self._size
                # print("from line " + str(row-wiskersLen) + "to line" + str(row+wiskersLen+1))
                # print(average)
                X[row, :] = average

                return X

    def _sizeSum(self, X, start, end):
        featureNumber = X.shape[1]
        result = np.zeros(featureNumber)
        for row in X[start:end, :]:
            result += row
        return result