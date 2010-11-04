""" Transformers to perform common preprocessing steps.
"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD

import numpy as np

from ..base import BaseEstimator


def _mean_and_std(X, axis=0, with_std=True):
    Xr = np.rollaxis(X, axis)
    mean_ = Xr.mean(axis=0)
    std_ = Xr.std(axis=0) if with_std else None
    return mean_, std_


def scale(X, axis=0, with_std=True, copy=True):
    """Method to standardize a dataset along any axis
    """
    mean_, std_ = _mean_and_std(X, axis, with_std)
    if copy:
        X = X.copy()
    Xr = np.rollaxis(X, axis)
    Xr -= mean_
    if with_std:
        std_[std_==0.0] = 1
        Xr /= std_
    return X


class Scaler(BaseEstimator):
    """Object to standardize a dataset

    It centers the dataset and optionaly scales to
    fix the variance to 1 for each feature

    """
    def __init__(self, with_std=True):
        self.with_std = with_std

    def fit(self, X, y=None, **params):
        self._set_params(**params)
        self.mean_, self.std_ = _mean_and_std(X, axis=0, with_std=self.with_std)
        return self

    def transform(self, X, y=None, copy=True):
        if copy:
            X = X.copy()
        # We are taking a view of the X array and modifying it
        X -= self.mean_
        if self.with_std:
            X /= self.std_
        return X

class Normalizer(BaseEstimator):
    """
    Normalize vectors such that they sum to 1.
    """

    def fit(self, X, y=None, **params):
        self._set_params(**params)
        return self

    def transform(self, X, y=None, copy=True):
        if copy:
            X = X.copy()

        X /= X.sum(axis=1)[:,np.newaxis]

        return X

class LengthNormalizer(BaseEstimator):
    """
    Normalize vectors to unit vectors.
    """

    def fit(self, X, y=None, **params):
        self._set_params(**params)
        return self

    def transform(self, X, y=None, copy=True):
        if copy:
            X = X.copy()

        X /= np.sqrt(np.sum(X ** 2, axis=1))[:,np.newaxis]

        return X

class Binarizer(BaseEstimator):
    """
    Binarize data according to a threshold.
    """

    def __init__(self, threshold=0.0):
        self.threshold = threshold

    def fit(self, X, y=None, **params):
        self._set_params(**params)
        return self

    def transform(self, X, y=None, copy=True):
        if copy:
            X = X.copy()

        cond = X > self.threshold
        not_cond = np.logical_not(cond)
        X[cond] = 1
        X[not_cond] = 0

        return X
