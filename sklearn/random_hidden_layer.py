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

from abc import ABCMeta, abstractmethod

import numpy as np
import scipy.sparse as sp

from math import sqrt
from scipy.spatial.distance import cdist

from .metrics import pairwise_distances
from .utils import check_random_state, atleast2d_or_csr
from .utils.extmath import safe_sparse_dot
from .base import BaseEstimator, TransformerMixin

__all__ = ['SimpleRandomHiddenLayer',
           'RBFRandomHiddenLayer']


# little function to propagate docstrings
# tinily tweaked from Paul McGuire
def _take_docstring_from(cls):
    def docstring_decorator(fn):
        fn.__doc__ = getattr(cls, fn.__name__).__doc__
        return fn
    return docstring_decorator


# Abstract Base Class for random hidden layers
class BaseRandomHiddenLayer(BaseEstimator, TransformerMixin):
    __metaclass__ = ABCMeta

    # take n_hidden and random_state, init components_ and
    # input_activations_
    def __init__(self, n_hidden=20, random_state=0, user_func=None,
                 user_args={}):

        self.n_hidden = n_hidden
        self.random_state = random_state
        self.user_func = user_func
        self.user_args = user_args

        self.components_ = dict()
        self.input_activations_ = None

    @abstractmethod
    def generate_components(self, X):
        """Generate components of hidden layer given X"""

    @abstractmethod
    def compute_input_activations(self, X):
        """Compute input activations given X"""

    @abstractmethod
    def compute_hidden_activations(self, X):
        """Compute hidden activations given X"""

    # perform fit by generating random components based
    # on the input array
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

        self.generate_components(X)

        return self

    # perform transformation by calling compute_hidden_activations
    # (which will normally call compute_input_activations first)
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

        if (self.components_ is None):
            raise ValueError('No components initialized')

        return self.compute_hidden_activations(X)


class SimpleRandomHiddenLayer(BaseRandomHiddenLayer):
    """Simple Random Hidden Layer transformer

    Creates a layer of units as a specified functions of an activation
    value determined by the dot product of the input and a random vector
    plus a random bias term:

     f(a), s.t. a = dot(x, hidden_weights) + bias

    and xfer function f() which defaults to numpy.tanh if not supplied,
    but can be any callable that returns an array of the same shape as
    its argument (input activation array, shape [n_samples, n_hidden])

    Parameters
    __________
    `n_hidden`      : int
        number of units to generate (default=20)

    `user_func` : callable
        function used to transform input activation (default=None)
        if not present or not callable, numpy.tanh is used instead

    `user_args` : dictionary
        contains keyword arguments for the transfer function

    `random_state`  : int, RandomState instance or None (default=None)
        Control the pseudo random number generator used to generate the
        hidden unit weights at fit time.

    Attributes
    __________
    `components_`  : dictionary containing two keys:
        `bias_weights_`   : numpy array of shape [n_hidden]
        `hidden_weights_` : numpy array of shape [n_features, n_hidden]

    `xfer_func_` : callable
        user_func if provided, called with keyword args from user_args
        dictionary, otherwise np.tanh

    See Also
    ________
    ELMRegressor, ELMClassifier, SimpleELMRegressor, SimpleELMClassifier,
    RBFRandomHiddenLayer
    """

    # default setup, plus initialization of xfer_func
    def __init__(self, n_hidden=20, random_state=None,
                 user_func=None, user_args={}):

        super(SimpleRandomHiddenLayer, self).__init__(n_hidden, random_state,
                                                      user_func, user_args)

    @_take_docstring_from(BaseRandomHiddenLayer)
    def generate_components(self, X):
        rand_state = check_random_state(self.random_state)
        n_features = X.shape[1]

        b_size = self.n_hidden
        hw_size = (n_features, self.n_hidden)

        self.components_['biases'] = rand_state.normal(size=b_size)
        self.components_['weights'] = rand_state.normal(size=hw_size)

    @_take_docstring_from(BaseRandomHiddenLayer)
    def compute_input_activations(self, X):
        b = self.components_['biases']
        w = self.components_['weights']

        self.input_activations_ = safe_sparse_dot(X, w)
        self.input_activations_ += b

    @_take_docstring_from(BaseRandomHiddenLayer)
    def compute_hidden_activations(self, X):
        self.compute_input_activations(X)

        if (not callable(self.user_func)):
            X_new = np.tanh(self.input_activations_)
        else:
            X_new = self.user_func(self.input_activations_, **self.user_args)

        return X_new


class RBFRandomHiddenLayer(BaseRandomHiddenLayer):
    """Random RBF Hidden Layer transformer

    Creates a layer of radial basis function units where:

       f(a), s.t. a = ||x-c||/r

    with c the unit center and r = max(||x-c||)/sqrt(n_centers*2).

    f() defaults to exp(-gamma * a^2)
    gamma defaults to 0.1

    If centers are not provided and use_exemplars is False (see below),
    then centers are uniformly distributed over the input space.

    Parameters
    __________
    `n_hidden`      : int
        optional (default=20)
        number of units to generate, ignored if centers are provided

    `centers`       : array of shape (n_hidden, n_features)
        optional (default=None)

        if provided, overrides various ways of internally computing the centers

    `radii`         : array of shape (n_hidden)
        optional (default=None)

        if provided, overrides internal radius calculation

    `use_exemplars` : bool
        optional (default=False)

        if True, uses random examples from the input to determine the RBF
        centers, ignored if centers are provided

    `use_radii` : bool
        optional (default=True)

        if True, divides the raw distances by the radii

    `gamma` : float
        optional (default=1.0)
        width multiplier for RBF argument

    `user_func` : callable
        function used to transform input activation (default=None)
        if not present or not callable, numpy.tanh is used instead

    `user_args` : dictionary
        contains keyword arguments for the transfer function

    `random_state`  : int or RandomState instance
        optional (default=None)

        Control the pseudo random number generator used to generate the
        centers at fit time, ignored if centers are provided

    Attributes
    __________
    `components_`   : dictionary containing two keys:
        `radii_`   : numpy array of shape [n_hidden]
        `centers_` : numpy array of shape [n_hidden, n_features]

    `xfer_func_` : callable
        user_func if provided, called with keyword args from user_args
        dictionary, otherwise

    See Also
    ________
    ELMRegressor, ELMClassifier, SimpleELMRegressor, SimpleELMClassifier,
    SimpleRandomHiddenLayer
    """

    def __init__(self, n_hidden=20, centers=None, radii=None, gamma=1.0,
                 use_exemplars=False, random_state=None, user_func=None,
                 user_args={}):

        super(RBFRandomHiddenLayer, self).__init__(n_hidden, random_state,
                                                   user_func, user_args)

        self.radii = radii
        self.centers = centers
        self.gamma = gamma
        self.use_exemplars = use_exemplars

    @_take_docstring_from(BaseRandomHiddenLayer)
    def generate_components(self, X):
        sparse = sp.issparse(X)
        self._compute_centers(X, sparse)
        self._compute_radii(X, sparse)

    @_take_docstring_from(BaseRandomHiddenLayer)
    def compute_input_activations(self, X):
        radii = self.components_['radii']
        centers = self.components_['centers']
        self.input_activations_ = cdist(X, centers)/radii

    @_take_docstring_from(BaseRandomHiddenLayer)
    def compute_hidden_activations(self, X):
        self.compute_input_activations(X)

        if (not callable(self.user_func)):
            X_new = np.exp(-self.gamma * pow(self.input_activations_, 2.0))
        else:
            X_new = self.user_func(self.input_activations_, **self.user_args)

        return X_new

    # determine centers
    def _compute_centers(self, X, sparse):
        # use given centers
        if (self.centers is not None):
            centers = self.centers

        else:
            n_samples, n_features = X.shape
            rs = check_random_state(self.random_state)

            # use examples from the data as centers
            if (self.use_exemplars):
                if (n_samples < self.n_hidden):
                    msg = "n_hidden must be <= n_samples when using exemplars"
                    raise ValueError(msg)
                max_index = X.shape[0] - 1
                indices = rs.permutation(max_index)[:self.n_hidden]
                centers = X[indices, :]

            # use uniformly distributed points from the input space as centers
            else:
                if (sparse):
                    min_X, max_X = min(X.data), max(X.data)
                else:
                    min_X, max_X = np.min(X), np.max(X)
                ctrs_size = (self.n_hidden, n_features)
                centers = rs.uniform(min_X, max_X, ctrs_size)

        self.components_['centers'] = centers

    # compute radii
    def _compute_radii(self, X, sparse):
        if (self.radii is not None):
            radii = self.radii

        else:
            centers = self.components_['centers']

            n_centers = centers.shape[0]
            max_dist = np.max(pairwise_distances(centers))
            radii = np.ones(n_centers) * max_dist/sqrt(2.0 * n_centers)

        self.components_['radii'] = radii
