# -*- coding: utf8
"""Random Hidden Layer transformers

Random hidden layers are arrays of hidden unit activations that are
random functions of input activation values (dot products for simple
activation functions, distances from prototypes for radial basis
functions).

They are used in the implementation of Extreme Learning Machines (ELMs),
but can be used as a general input mapping.
"""

# Author: David C. Lambert <dcl@panix.com>
# License: Simple BSD

import numpy as np
import scipy.sparse as sp

from math import sqrt
from scipy.spatial.distance import cdist

from .cluster import k_means
from .metrics import pairwise_distances
from .utils import check_random_state, atleast2d_or_csr
from .utils.extmath import safe_sparse_dot
from .base import BaseEstimator, TransformerMixin

__all__ = ['SimpleRandomHiddenLayer',
           'RBFRandomHiddenLayer']


class SimpleRandomHiddenLayer(BaseEstimator, TransformerMixin):
    """Random Hidden Layer transformer

    Creates a layer of units as a specified functions of an activation
    value determined by the dot product of the input and a random vector
    plus a random bias term:

     f(a), s.t. a = dot(x, hidden_weights) + bias

    and transfer function f() which defaults to numpy.tanh if not supplied,
    but can be any callable that returns an array of the same shape as
    its argument (input activation array, shape [n_samples, n_hidden])


    Parameters
    __________
    `n_hidden`      : int
        number of units to generate (default=20)

    `transfer_func` : callable
        function used to transform input activation (default=None)
        if not present or not callable, numpy.tanh is used instead

    `random_state`  : int, RandomState instance or None (default=None)
        Control the pseudo random number generator used to generate the
        hidden unit weights at fit time.

    Attributes
    __________
    `bias_weights_`   : numpy array of shape [n_hidden]
    `hidden_weights_` : numpy array of shape [n_features, n_hidden]

    See Also
    ________
    ELMRegressor, ELMClassifier, SimpleELMRegressor, SimpleELMClassifier,
    RBFRandomHiddenLayer
    """

    def __init__(self, n_hidden=20, transfer_func=None, random_state=None):
        self.bias_weights_ = None
        self.hidden_weights_ = None

        self.n_hidden = n_hidden
        self.transfer_func = transfer_func
        self.random_state = random_state

    def fit(self, X, y=None):
        """Generate a random hidden layer.

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            Training set: only the shape is used to generate random hidden
            weights for random functional projections

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        self

        """
        X = atleast2d_or_csr(X)

        n_features = X.shape[1]

        rs = check_random_state(self.random_state)

        self.bias_weights_ = rs.normal(0.0, 1.0, self.n_hidden)

        hw_size = (n_features, self.n_hidden)
        self.hidden_weights_ = rs.normal(0.0, 1.0, hw_size)

        return self

    def transform(self, X, y=None):
        """Generate the random hidden layer's activations given X as input.

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            The input data to project into a random hidden layer's activations

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        X_new : numpy array or scipy sparse of shape [n_samples, n_components]
            Projected array.

        """
        X = atleast2d_or_csr(X)

        if (self.hidden_weights_ is None):
            raise ValueError('No random weights initialized')

        input_activation = safe_sparse_dot(X, self.hidden_weights_)
        input_activation += self.bias_weights_

        if (not callable(self.transfer_func)):
            X_new = np.tanh(input_activation)
        else:
            X_new = self.transfer_func(input_activation)

        return X_new


class RBFRandomHiddenLayer(BaseEstimator, TransformerMixin):
    """Random RBF Hidden Layer transformer

    Creates a layer of radial basis function units where:

       f(x) = exp(-gamma * (||x-c||/r)^2)

    with c the unit center and r = max(||x-c||)/sqrt(n_centers*2).

    If centers are not provided and neither use_kmeans nor use_exemplars
    is True (see below), centers are uniformly distributed over the input
    space.

    Parameters
    __________
    `n_hidden`      : int
        optional (default=20)
        number of units to generate, ignored if centers are provided

    `gamma`         : float
        optional (default=1.0)
        width multiplier for RBF argument

    `centers`       : array of shape (n_hidden, n_features)
        optional(default=None)

        if provided, overrides various ways of internally computing the centers

    `use_kmeans`    : bool
        optional (default=False)

        if True, uses vanilla k_means from the input data to determine the RBF
        centers, ignored if centers are provided

    `use_exemplars` : bool
        optional (default=False)

        if True, uses random examples from the input to determine the RBF
        centers, ignored if centers are provided

    `random_state`  : int or RandomState instance
        optional (default=None)

        Control the pseudo random number generator used to generate the
        centers at fit time, ignored if centers are provided

    Attributes
    __________
    `unit_radii_`   : numpy array of shape [n_hidden]
    `unit_centers_` : numpy array of shape [n_hidden, n_features]
    `centers`       : numpy array of shape [n_hidden, n_features]
    `use_kmeans`    : bool, flag set when using kmeans
    `use_exemplars` : bool, flag set when using exemplars

    See Also
    ________
    ELMRegressor, ELMClassifier, SimpleELMRegressor, SimpleELMClassifier,
    SimpleRandomHiddenLayer
    """

    def __init__(self, n_hidden=20, gamma=1.0, centers=None,
                 use_kmeans=False, use_exemplars=False, random_state=None):

        self.unit_radii_ = None
        self.unit_centers_ = None

        self.n_hidden = n_hidden
        self.gamma = gamma
        self.centers = centers
        self.use_kmeans = use_kmeans
        self.use_exemplars = use_exemplars
        self.random_state = random_state

    # determine centers
    def _get_centers(self, X, sparse):
        # centers given to __init__ already
        if (self.centers is not None):
            self.unit_centers_ = self.centers
            return

        n_samples, n_features = X.shape

        rs = check_random_state(self.random_state)

        # use kmeans to determine centers
        if (self.use_kmeans):
            if (n_samples < self.n_hidden):
                msg = "n_hidden must be <= n_samples when using kmeans"
                raise ValueError(msg)
            (self.unit_centers_, _, _) = k_means(X, self.n_hidden)

        # use examples from the data as centers
        elif (self.use_exemplars):
            if (n_samples < self.n_hidden):
                msg = "n_hidden must be <= n_samples when using exemplars"
                raise ValueError(msg)
            max_index = X.shape[0] - 1
            indices = rs.permutation(max_index)[:self.n_hidden]
            self.unit_centers_ = X[indices, :]

        # use uniformly distributed points from the input space as centers
        else:
            if (sparse):
                min_X, max_X = min(X.data), max(X.data)
            else:
                min_X, max_X = np.min(X), np.max(X)
            ctrs_size = (self.n_hidden, n_features)
            self.unit_centers_ = rs.uniform(min_X, max_X, ctrs_size)

    # compute radii
    def _compute_radii(self, X, sparse):
        n_centers = self.unit_centers_.shape[0]
        max_dist = np.max(pairwise_distances(self.unit_centers_))
        self.unit_radii_ = np.ones(n_centers) * max_dist/sqrt(2.0*n_centers)

    # initialize centers and radii given X
    def fit(self, X, y=None):
        """Generate a random hidden layer.

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            Training set: only the shape is used to generate random hidden
            weights for random functional projections

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        self

        """
        X = atleast2d_or_csr(X)

        sparse = sp.issparse(X)
        self._get_centers(X, sparse)
        self._compute_radii(X, sparse)
        return self

    # compute RBF activations given X
    def transform(self, X, y=None):
        """Generate the random hidden layer's activations given X as input.

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            The input data to project into a random hidden layer's activations

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        X_new : numpy array or scipy sparse of shape [n_samples, n_components]
            Projected array.

        """
        X = atleast2d_or_csr(X)

        distances = cdist(X, self.unit_centers_)
        X_new = np.exp(-self.gamma * pow(distances/self.unit_radii_, 2.0))
        return X_new
