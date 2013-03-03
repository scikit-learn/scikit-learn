# -*- coding: utf8
"""The :mod:`sklearn.neural_networks.random_hidden_layer` module
implements Random Hidden Layer transformers.

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

from math import sqrt

import numpy as np
import scipy.sparse as sp
from scipy.spatial.distance import cdist

from ..metrics import pairwise_distances
from ..utils import check_random_state, atleast2d_or_csr
from ..utils.extmath import safe_sparse_dot
from ..base import BaseEstimator, TransformerMixin

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

    _internal_xfer_funcs = dict()

    # take n_hidden and random_state, init components_ and
    # input_activations_
    def __init__(self, n_hidden=20, random_state=0, xfer_func=None,
                 xfer_args=None):

        self.n_hidden = n_hidden
        self.random_state = random_state
        self.xfer_func = xfer_func
        self.xfer_args = xfer_args

        self.components_ = dict()
        self.input_activations_ = None

        # keyword args for internally defined funcs
        self._extra_args = dict()

    @abstractmethod
    def _generate_components(self, X):
        """Generate components of hidden layer given X"""

    @abstractmethod
    def _compute_input_activations(self, X):
        """Compute input activations given X"""

    # compute input activations and pass them
    # through the hidden layer transfer functions
    # to compute the transform
    def _compute_hidden_activations(self, X):
        """Compute hidden activations given X"""

        self._compute_input_activations(X)

        acts = self.input_activations_

        if (callable(self.xfer_func)):
            args_dict = self.xfer_args if (self.xfer_args) else dict()
            X_new = self.xfer_func(acts, **args_dict)
        else:
            func_name = self.xfer_func
            func = self._internal_xfer_funcs[func_name]

            X_new = func(acts, **self._extra_args)

        return X_new

    # perform fit by generating random components based
    # on the input array
    def fit(self, X, y=None):
        """Generate a random hidden layer.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape [n_samples, n_features]
            Training set: only the shape is used to generate random component
            values for hidden units

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        self
        """
        X = atleast2d_or_csr(X)

        self._generate_components(X)

        return self

    # perform transformation by calling compute_hidden_activations
    # (which will normally call compute_input_activations first)
    def transform(self, X, y=None):
        """Generate the random hidden layer's activations given X as input.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape [n_samples, n_features]
            Data to transform

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        X_new : numpy array of shape [n_samples, n_components]
        """
        X = atleast2d_or_csr(X)

        if (self.components_ is None):
            raise ValueError('No components initialized')

        return self._compute_hidden_activations(X)


class SimpleRandomHiddenLayer(BaseRandomHiddenLayer):
    """Simple Random Hidden Layer transformer

    Creates a layer of units as a specified functions of an activation
    value determined by the dot product of the input and a random vector
    plus a random bias term:

     f(a), s.t. a = dot(x, hidden_weights) + bias

    and transfer function f() which defaults to numpy.tanh if not supplied
    but can be any callable that returns an array of the same shape as
    its argument (input activation array, shape [n_samples, n_hidden])

    Parameters
    ----------
    `n_hidden` : int, optional (default=20)
        Number of units to generate

    `xfer_func` : {callable, string} optional (default='tanh')
        Function used to transform input activation
        It must be one of 'tanh', 'sine', 'tribas', 'sigmoid', 'hardlim' or
        a callable.  If none is given, 'tanh' will be used. If a callable
        is given, it will be used to compute the hidden unit activations.

    `xfer_args` : dictionary, optional (default=None)
        Supplies keyword arguments for a callable xfer_func

    `random_state`  : int, RandomState instance or None (default=None)
        Control the pseudo random number generator used to generate the
        hidden unit weights at fit time.

    Attributes
    ----------
    `input_activations_` : numpy array of shape [n_samples, n_hidden]
        Array containing dot(x, hidden_weights) + bias for all samples

    `components_` : dictionary containing two keys:
        `bias_weights_`   : numpy array of shape [n_hidden]
        `hidden_weights_` : numpy array of shape [n_features, n_hidden]

    See Also
    --------
    ELMRegressor, ELMClassifier, SimpleELMRegressor, SimpleELMClassifier,
    RBFRandomHiddenLayer
    """

    #
    # internal transfer function (RBF) definitions
    #

    # triangular transfer function
    _tribas = (lambda x: np.clip(1.0 - np.fabs(x), 0.0, 1.0))

    # sigmoid transfer function
    _sigmoid = (lambda x: 1.0/(1.0 + np.exp(-x)))

    # hard limit transfer function
    _hardlim = (lambda x: np.array(x > 0.0, dtype=float))

    # internal transfer function table
    _internal_xfer_funcs = {'sine': np.sin,
                            'tanh': np.tanh,
                            'tribas': _tribas,
                            'sigmoid': _sigmoid,
                            'hardlim': _hardlim
                            }

    # default setup, plus initialization of xfer_func
    def __init__(self, n_hidden=20, random_state=None,
                 xfer_func='tanh', xfer_args=None):

        super(SimpleRandomHiddenLayer, self).__init__(n_hidden, random_state,
                                                      xfer_func, xfer_args)

        if (isinstance(self.xfer_func, str)):
            if (self.xfer_func not in self._internal_xfer_funcs.keys()):
                msg = "unknown transfer function '{}'".format(self.xfer_func)
                raise ValueError(msg)

    @_take_docstring_from(BaseRandomHiddenLayer)
    def _generate_components(self, X):
        rand_state = check_random_state(self.random_state)
        n_features = X.shape[1]

        b_size = self.n_hidden
        hw_size = (n_features, self.n_hidden)

        self.components_['biases'] = rand_state.normal(size=b_size)
        self.components_['weights'] = rand_state.normal(size=hw_size)

    @_take_docstring_from(BaseRandomHiddenLayer)
    def _compute_input_activations(self, X):
        b = self.components_['biases']
        w = self.components_['weights']

        self.input_activations_ = safe_sparse_dot(X, w)
        self.input_activations_ += b


# Random Hidden Layer of radial basis function units
class RBFRandomHiddenLayer(BaseRandomHiddenLayer):
    """Random RBF Hidden Layer transformer

    Creates a layer of radial basis function units where:

       f(a), s.t. a = ||x-c||/r

    with c the unit center and r = max(||x-c||)/sqrt(n_centers*2).

    f() defaults to exp(-gamma * a^2) (gaussian rbf)
    gamma defaults to 1.0

    If centers are not provided and use_exemplars is False (see below),
    then centers are uniformly distributed over the input space.

    Parameters
    ----------
    `n_hidden` : int, optional (default=20)
        Number of units to generate, ignored if centers are provided

    `xfer_func` : {callable, string} optional (default='gaussian')
        Function used to transform input activation.
        It must be one of 'gaussian', 'poly_spline', 'multiquadric' or
        a callable.  If none is given, 'gaussian' will be used. If a
        callable is given, it will be used to compute the hidden unit
        activations.

    `xfer_args` : dictionary, optional (default=None)
        Supplies keyword arguments for a callable xfer_func

    `gamma` : {int, float} optional (default=1.0)
        Width multiplier for RBF distance argument, ignored if callable
        xfer_func is provided.  Must be an int > 0 when xfer_func is
        'poly_spline'.

    `centers` : array of shape (n_hidden, n_features), optional (default=None)
        If provided, overrides internal computation of the centers

    `radii` : array of shape (n_hidden),  optional (default=None)
        If provided, overrides internal computation of the radii

    `use_exemplars` : bool, optional (default=False)
        If True, uses random examples from the input to determine the RBF
        centers, ignored if centers are provided

    `random_state`  : int or RandomState instance, optional (default=None)
        Control the pseudo random number generator used to generate the
        centers at fit time, ignored if centers are provided

    Attributes
    ----------
    `components_` : dictionary containing two keys:
        `radii_`   : numpy array of shape [n_hidden]
        `centers_` : numpy array of shape [n_hidden, n_features]

    `input_activations_` : numpy array of shape [n_samples, n_hidden]
        Array containing ||x-c||/r for all samples

    See Also
    --------
    ELMRegressor, ELMClassifier, SimpleELMRegressor, SimpleELMClassifier,
    SimpleRandomHiddenLayer
    """

    #
    # internal transfer function (RBF) definitions
    #

    # gaussian RBF
    _gaussian = (lambda x, gamma: np.exp(-gamma * pow(x, 2.0)))

    # multiquadric spline RBF
    _multiquadric = (lambda x, gamma:
                         np.sqrt(1.0 + pow(gamma * x, 2.0)))

    # polyharmonic spline RBF
    def _poly_spline(acts, gamma):
        if (not isinstance(gamma, int) or gamma < 1):
            msg = 'Gamma must be integer > 0 for poly_spline'
            raise ValueError(msg)

        # add epsilon to avoid log(0) exception
        epsilon = 1.0e-8
        acts += epsilon

        X_new = pow(acts, gamma)
        if ((gamma % 2) == 0):
            X_new *= np.log(acts)

        return X_new

    # internal RBF table
    _internal_xfer_funcs = {'gaussian': _gaussian,
                            'poly_spline': _poly_spline,
                            'multiquadric': _multiquadric
                            }

    def __init__(self, n_hidden=20, random_state=None, xfer_func='gaussian',
                 xfer_args=None, gamma=1.0, centers=None, radii=None,
                 use_exemplars=False):

        super(RBFRandomHiddenLayer, self).__init__(n_hidden, random_state,
                                                   xfer_func, xfer_args)

        if (isinstance(self.xfer_func, str)):
            if (self.xfer_func not in self._internal_xfer_funcs.keys()):
                msg = "unknown transfer function '{}'".format(self.xfer_func)
                raise ValueError(msg)

        self.radii = radii
        self.centers = centers
        self.gamma = gamma
        self.use_exemplars = use_exemplars

    # property methods for 'gamma' arg, use
    # self._extra_args dictionary
    @property
    def gamma(self):
        return self._extra_args['gamma']

    @gamma.setter
    def gamma(self, value):
        self._extra_args['gamma'] = value

    @_take_docstring_from(BaseRandomHiddenLayer)
    def _generate_components(self, X):
        sparse = sp.issparse(X)
        self._compute_centers(X, sparse)
        self._compute_radii(X, sparse)

    @_take_docstring_from(BaseRandomHiddenLayer)
    def _compute_input_activations(self, X):
        radii = self.components_['radii']
        centers = self.components_['centers']

        self.input_activations_ = cdist(X, centers)/radii

    # determine centers
    def _compute_centers(self, X, sparse):
        # use supplied centers
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

                max_index = n_samples - 1
                indices = rs.permutation(max_index)[:self.n_hidden]
                centers = X[indices, :]

            # use uniformly distributed points from the input space as centers
            else:
                if (sparse):
                    X_dtype = X.dtype.type(0)
                    min_X = np.minimum(X_dtype, np.min(X.data))
                    max_X = np.maximum(X_dtype, np.max(X.data))
                else:
                    min_X, max_X = np.min(X), np.max(X)

                ctrs_size = (self.n_hidden, n_features)
                centers = rs.uniform(min_X, max_X, ctrs_size)

        self.components_['centers'] = centers

    # compute radii
    def _compute_radii(self, X, sparse):
        # use supplied radii
        if (self.radii is not None):
            radii = self.radii

        else:
            centers = self.components_['centers']

            n_centers = centers.shape[0]
            max_dist = np.max(pairwise_distances(centers))
            radii = np.ones(n_centers) * max_dist/sqrt(2.0 * n_centers)

        self.components_['radii'] = radii
