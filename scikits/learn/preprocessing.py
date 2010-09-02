import numpy as np

from .base import BaseEstimator

class Scaler(BaseEstimator):
    """Object to standardize a dataset along any axis
    
    It centers the dataset and optionaly scales to 
    fix the variance to 1.
    
    """
    def __init__(self, axis=0, with_std=True):
        self.axis = axis
        self.with_std = with_std

    def fit(self, X, y=None, **params):
        self._set_params(**params)
        X = np.rollaxis(X, self.axis)
        self.mean = X.mean(axis=0)
        if self.with_std:
            self.std = X.std(axis=0)
        return self

    def transform(self, X, y=None, copy=True):
        if copy:
            X = X.copy()
        Xr = np.rollaxis(X, self.axis)
        Xr -= self.mean
        if self.with_std:
            Xr /= self.std
        return X
