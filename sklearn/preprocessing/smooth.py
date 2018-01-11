# Authors: Luigi De Bianchi <luigi.debianchi@gmail.com>
# License: BSD 3 clause

from __future__ import division

import numpy as np
from scipy import sparse

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.validation import (check_is_fitted, FLOAT_DTYPES)


__all__ = [
    'RectangularSmoother',
]


class RectangularSmoother(BaseEstimator, TransformerMixin):

    _modes = [ 'mirror', 'constant', 'nearest', 'wrap', 'interp']
    # TODO: implements interp

    def __init__(self, size=3, mode='constant', cval=0.0):
        if mode not in self._modes:
            raise ValueError("Given mode is not supported.")

        if size%2 == 0:
            raise ValueError("Size must be odd.")

        self._mode = mode
        self._size = size
        self._whisker = int((self._size - 1) / 2)
        self._cval = cval

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

        array_filler_up, array_filler_down = self._populate_filler((self._whisker, X.shape[1]))

        supported_input = np.concatenate((array_filler_up, X, array_filler_down), axis=0)

        result = np.zeros(X.shape)
        result[0, :] = self.sum_lines(0, self._size, )
        for row in range(1, X.shape[0]):
            result[row, :] = result[row - 1, :] - supported_input[row - 1, :] + supported_input[row + self._size - 1, :]
        result = np.divide(result, self._size)
        return result

    @staticmethod
    def sum_lines(X, start, end):
        result = np.zeros(X.shape[1])
        for row in X[start:end, :]:
            result += row
        return result

    def _populate_filler(self, num, X):
        filler_up = np.zeros((num, X.shape[1]))
        filler_down = np.zeros((num, X.shape[1]))

        if self._mode == 'mirror':
            return filler_up, filler_down

        if self._mode == 'constant':
            filler_up[:, :] = self._cval
            filler_down[:, :] = self._cval
            return filler_up, filler_down

        if self._mode == 'nearest':
            filler_up[:, :] = X[0, :]
            filler_down[:, :] = X[-1, :]
            return filler_up, filler_down

        if self._mode == 'wrap':
            return filler_up, filler_down

