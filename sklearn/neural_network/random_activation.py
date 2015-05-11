from abc import ABCMeta, abstractmethod

import numpy as np

from sklearn.base import TransformerMixin
from sklearn.externals import six
from sklearn.utils import check_random_state
from sklearn.utils.extmath import safe_sparse_dot
from sklearn.utils.validation import check_array


class RandomActivation(six.with_metaclass(ABCMeta, TransformerMixin)):
    def __init__(self, n_activated_features=10, weight_scale='auto', 
                 activation=False, intercept=True, random_state=None):
        self.n_activated_features = n_activated_features
        self.weight_scale = weight_scale
        self.intercept = intercept
        self.random_state = random_state

    def fit(self, X, y=None):
        """
        Generate random parameters based on the number of features the input
        data has.

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features].

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        self

        """
        X = check_array(X, accept_sparse=['csr', 'csc'])

        n_samples, n_features = X.shape
        
        rng = check_random_state(self.random_state)

        if self.weight_scale == 'auto':
            weight_init_bound = np.sqrt(6. / (n_features + 
                                              self.n_activated_features))
        else:
            weight_init_bound = self.weight_scale

        self.coef_ = rng.uniform(-weight_init_bound, weight_init_bound,
                                (n_features, self.n_activated_features))

        if self.intercept:
            self.intercept_ = rng.uniform(-weight_init_bound, weight_init_bound,
                                           self.n_activated_features)

        return self

    def transform(self, X, y=None):
        """
        Transform input data to another space using the randomly generated 
        parameters

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features].

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        X_new : numpy array or scipy sparse of shape [n_samples, n_components]
            Projected array.

        """
        X = check_array(X, accept_sparse=['csr', 'csc'])

        X_new = safe_sparse_dot(X, self.coef_)

        if self.intercept:
            X_new += self.intercept_

        return X_new