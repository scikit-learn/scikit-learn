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

    _modes = ['mirror', 'constant', 'nearest', 'wrap', 'interp']

    def __init__(self, size=3, mode='constant', cval=0.0):
        if mode not in self._modes:
            raise ValueError("Given mode is not supported.")

        if size % 2 == 0:
            raise ValueError("Size must be odd.")

        self.mode = mode
        self.size = size
        self.cval = cval

    def _reset(self):
        if hasattr(self, '_mode'):
            del self.mode
            del self.size
            del self._whisker
            del self.cval

    def fit(self, x, y=None):
        if sparse.issparse(x):
            raise TypeError("Smoother does no support sparse input. ")

        x = check_array(x, copy=True, warn_on_dtype=True,
                        estimator=self, dtype=FLOAT_DTYPES)

        self._whisker = int((self.size - 1) / 2)

        return self

    def transform(self, x):
        check_is_fitted(self, '_whisker')

        x = check_array(x, copy=True, dtype=FLOAT_DTYPES)

        array_filler_up, array_filler_down = self._populate_fillers(x)

        supported_input = np.concatenate((array_filler_up, x, array_filler_down), axis=0)
        result = np.zeros(x.shape)

        if self.mode == 'interp':
            supported_input = x
            result = np.zeros((x.shape[0] - self._whisker, x.shape[1]))

        result[0, :] = self.sum_lines(supported_input, 0, self.size)
        for row in range(1, result.shape[0]):
            result[row, :] = result[row - 1, :] - supported_input[row - 1, :] + supported_input[row + self.size - 1, :]
        result = np.divide(result, self.size)
        return result

    @staticmethod
    def sum_lines(x, start, end):
        result = np.zeros(x.shape[1])
        for row in x[start:end, :]:
            result += row
        return result

    def _populate_fillers(self, x):

        if x.shape[0] < self._whisker:
            raise TypeError("Too few sample with respect to the chosen window size")

        filler_up = np.zeros((self._whisker, x.shape[1]))
        filler_down = np.zeros((self._whisker, x.shape[1]))

        if self.mode == 'mirror':
            for i in range(0, self._whisker):
                filler_up[i, :] = x[self._whisker - i, :]
                filler_down[i, :] = x[- 2 - i, :]
            return filler_up, filler_down

        if self.mode == 'constant':
            filler_up[:, :] = self.cval
            filler_down[:, :] = self.cval
            return filler_up, filler_down

        if self.mode == 'nearest':
            filler_up[:, :] = x[0, :]
            filler_down[:, :] = x[-1, :]
            return filler_up, filler_down

        if self.mode == 'wrap':
            filler_up[:, :] = x[-self._whisker:, :]
            filler_down[:, :] = x[:self._whisker, :]
            return filler_up, filler_down
