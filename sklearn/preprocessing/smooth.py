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

    def fit(self, x, y=None):
        if sparse.issparse(x):
            raise TypeError("Smoother does no support sparse input. ")

        x = check_array(x, ensure_min_samples=int((self.size - 1) / 2), warn_on_dtype=True, estimator=self)
        self._feature_indices = x.shape[1]
        return self

    def transform(self, x):
        check_is_fitted(self, '_whisker')

        whisker = int((self.size - 1) / 2)
        x = check_array(x, ensure_min_samples=whisker, warn_on_dtype=True, estimator=self)

        if x.shape[1] != self._feature_indices:
            raise ValueError("X has different shape than during fitting."
                             " Expected %d, got %d."
                             % (self._feature_indices, x.shape[1]))

        array_filler_up, array_filler_down = self.populate_fillers(x, self.mode, whisker, self.cval)

        supported_input = np.concatenate((array_filler_up, x, array_filler_down), axis=0)

        if self.mode == 'interp':
            result = np.zeros((x.shape[0] - whisker, x.shape[1]))
        else:
            result = np.zeros(x.shape)

        result[0, :] = self.sum_samples(supported_input, 0, self.size)
        for row in range(1, result.shape[0]):
            result[row, :] = result[row - 1, :] - supported_input[row - 1, :] + supported_input[row + self.size - 1, :]
        result = np.divide(result, self.size)
        return result

    @staticmethod
    def sum_samples(x, start, end):
        result = np.zeros(x.shape[1])
        for row in x[start:end, :]:
            result += row
        return result

    @staticmethod
    def populate_fillers(x, mode, whisker, cval=None):

        if x.shape[0] < whisker:
            raise ValueError("Too few sample with respect to the chosen window size")

        if mode == 'interp':
            filler = np.zeros((0, x.shape[1]))
            return filler, filler

        filler_up = np.zeros((whisker, x.shape[1]))
        filler_down = np.zeros((whisker, x.shape[1]))

        if mode == 'mirror':
            for i in range(0, whisker):
                filler_up[i, :] = x[whisker - i, :]
                filler_down[i, :] = x[- 2 - i, :]
            return filler_up, filler_down

        if mode == 'constant':
            filler_up[:, :] = cval
            return filler_up, filler_up

        if mode == 'nearest':
            filler_up[:, :] = x[0, :]
            filler_down[:, :] = x[-1, :]
            return filler_up, filler_down

        if mode == 'wrap':
            filler_up[:, :] = x[-whisker:, :]
            filler_down[:, :] = x[:whisker, :]
            return filler_up, filler_down

