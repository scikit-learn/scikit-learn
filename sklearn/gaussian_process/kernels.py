""" Kernels for Gaussian process regression and classification.

The kernels in this module allow kernel-engineering, i.e., they can be
combined via the "+" and "*" operators. These expressions can also contain
scalar values, which are automatically converted to a constant kernel.

All kernels allow (analytic) gradient-based hyperparameter optimization.
The space of hyperparameters can be specified by giving lower und upper
boundaries for the value of each hyperparameter (the search space is thus
rectangular). This can be achieved by using a pair or triple instead of a
single float wherever a parameter value is specified. In case of a pair,
the first value specifies the lower boundary and the second value the upper
boundary. In case of a triple, the middle value specified the initial value
of the parameter during hyperparameter-optimization.
"""

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

# Note: this module is strongly inspired by the kernel module of the george
#       package.

from abc import ABCMeta, abstractmethod
from functools import partial

import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform
from scipy.optimize import approx_fprime

from ..metrics.pairwise import pairwise_kernels
from ..externals import six


class Kernel(six.with_metaclass(ABCMeta)):
    """ Base class for all kernels."""

    def _parse_param_space(self, param_space):
        if not np.iterable(param_space):  # fixed hyperparameter
            self.params = np.array([float(param_space)])
            self.has_bounds = False
            return
        param_space = np.atleast_2d(param_space)
        if param_space.shape[1] == 1:
            # fixed hyperparameter
            self.params = param_space[:, 0]
            self.has_bounds = False
        elif param_space.shape[1] == 2:
            # lower+upper bound for hyperparameter
            self.bounds = param_space
            self.has_bounds = True
            # Use geometric mean of upper and lower boundary as initial
            # hyperparameter value
            if np.any(np.equal(self.l_bound, None)) \
               or np.any(np.equal(self.u_bound, None)):
                raise ValueError("Lower or upper bound being None requires "
                                 "explicitly specifying the initial value.")
            self.params = np.array([np.sqrt(self.l_bound * self.u_bound)])
        elif param_space.shape[1] == 3:
            # lower bound, initial value, upper bound
            self.params = param_space[:, 1]
            self.bounds = param_space[:, [0, 2]]
            self.has_bounds = True
        else:
            raise ValueError("Invalid parameter space given. Must not have "
                             "more than 3 entries per parameter.")

    @property
    def n_params(self):
        return self.params.shape[0]

    @property
    def bounds(self):
        return np.vstack((self.l_bound, self.u_bound)).T

    @bounds.setter
    def bounds(self, bounds):
        bounds = bounds.reshape(-1, 2)
        self.l_bound = bounds[:, 0]
        self.u_bound = bounds[:, 1]

    def __add__(self, b):
        if not isinstance(b, Kernel):
            return Sum(self, ConstantKernel(b))
        return Sum(self, b)

    def __radd__(self, b):
        if not isinstance(b, Kernel):
            return Sum(ConstantKernel(b), self)
        return Sum(b, self)

    def __mul__(self, b):
        if not isinstance(b, Kernel):
            return Product(self, ConstantKernel(b))
        return Product(self, b)

    def __rmul__(self, b):
        if not isinstance(b, Kernel):
            return Product(ConstantKernel(b), self)
        return Product(b, self)

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__,
                                 ", ".join(map("{0}".format, self.params)))

    @abstractmethod
    def __call__(self, X, Y=None, eval_gradient=False):
        """Evaluate the kernel."""


class KernelOperator(Kernel):
    """ Base class for all kernel operators. """

    def __init__(self, k1, k2):
        self.k1 = k1
        self.k2 = k2
        # XXX: Deal with situations in which only some of the hyperparameter
        #      shall be optimized
        self.has_bounds = k1.has_bounds and k2.has_bounds

    @property
    def params(self):
        return np.append(self.k1.params, self.k2.params)

    @params.setter
    def params(self, theta):
        i = self.k1.n_params
        self.k1.params = theta[:i]
        self.k2.params = theta[i:]

    @property
    def bounds(self):
        assert self.has_bounds
        return np.vstack((self.k1.bounds, self.k2.bounds))

    @bounds.setter
    def bounds(self, bounds):
        i = self.k1.n_params
        self.k1.bounds = bounds[:i]
        self.k2.bounds = bounds[i:]


class Sum(KernelOperator):
    """ Sum-kernel k1 + k2 of two kernels k1 and k2.

    The resulting kernel is defined as
    k_sum(X, Y) = k1(X, Y) + k2(X, Y)

    Parameters
    ----------
    k1 : Kernel object
        The first base-kernel of the sum-kernel

    k2 : Kernel object
        The second base-kernel of the sum-kernel
    """

    def __call__(self, X, Y=None, eval_gradient=False):
        """ Return the kernel k(X, Y) and optionally its gradient.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Y : array, shape (n_samples_Y, n_features), (optional, default=None)
            Right argument of the returned kernel k(X, Y). If None, k(X, X)
            if evaluated instead.

        eval_gradient : bool (optional, default=False)
            Determines whether the gradient with respect to the kernel
            hyperparameter is determined.

        Returns
        -------
        K : array, shape (n_samples_X, n_samples_Y)
            Kernel k(X, Y)

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_params)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        if eval_gradient:
            K1, K1_gradient = self.k1(X, Y, eval_gradient=True)
            K2, K2_gradient = self.k2(X, Y, eval_gradient=True)
            return K1 + K2, np.dstack((K1_gradient, K2_gradient))
        else:
            return self.k1(X, Y) + self.k2(X, Y)

    def __repr__(self):
        return "{0} + {1}".format(self.k1, self.k2)


class Product(KernelOperator):
    """ Product-kernel k1 * k2 of two kernels k1 and k2.

    The resulting kernel is defined as
    k_prod(X, Y) = k1(X, Y) * k2(X, Y)

    Parameters
    ----------
    k1 : Kernel object
        The first base-kernel of the product-kernel

    k2 : Kernel object
        The second base-kernel of the product-kernel
    """

    def __call__(self, X, Y=None, eval_gradient=False):
        """ Return the kernel k(X, Y) and optionally its gradient.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Y : array, shape (n_samples_Y, n_features), (optional, default=None)
            Right argument of the returned kernel k(X, Y). If None, k(X, X)
            if evaluated instead.

        eval_gradient : bool (optional, default=False)
            Determines whether the gradient with respect to the kernel
            hyperparameter is determined.

        Returns
        -------
        K : array, shape (n_samples_X, n_samples_Y)
            Kernel k(X, Y)

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_params)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        if eval_gradient:
            K1, K1_gradient = self.k1(X, Y, eval_gradient=True)
            K2, K2_gradient = self.k2(X, Y, eval_gradient=True)
            return K1 * K2, np.dstack((K1_gradient * K2[:, :, np.newaxis],
                                       K2_gradient * K1[:, :, np.newaxis]))
        else:
            return self.k1(X, Y) * self.k2(X, Y)

    def __repr__(self):
        return "{0} * {1}".format(self.k1, self.k2)


class ConstantKernel(Kernel):
    """ Constant kernel.

    Can be used as part of a product-kernel where it scales the magnitude of
    the other factor (kernel) or as part of a sum-kernel, where it modifies
    the mean of the Gaussian process.
    """

    def __init__(self, param_space=1.0):
        self._parse_param_space(param_space)

    @property
    def params(self):
        return np.array([self.value])

    @params.setter
    def params(self, theta):
        assert len(theta) == 1
        self.value = theta[0]

    def __call__(self, X, Y=None, eval_gradient=False):
        """ Return the kernel k(X, Y) and optionally its gradient.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Y : array, shape (n_samples_Y, n_features), (optional, default=None)
            Right argument of the returned kernel k(X, Y). If None, k(X, X)
            if evaluated instead.

        eval_gradient : bool (optional, default=False)
            Determines whether the gradient with respect to the kernel
            hyperparameter is determined. Only supported when Y is None.

        Returns
        -------
        K : array, shape (n_samples_X, n_samples_Y)
            Kernel k(X, Y)

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_params)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        if Y is None:
            Y = X
        elif eval_gradient:
            raise ValueError("Gradient can only be evaluated when Y is None.")

        K = self.value * np.ones((X.shape[0], Y.shape[0]))
        if eval_gradient:
            return K, np.ones((X.shape[0], X.shape[0], 1))
        else:
            return K

    def __repr__(self):
        return "{0}".format(self.value)


class RBF(Kernel):
    """ Radial-basis function kernel (aka squared-exponential kernel).

    Both isotropic and anisotropic version are supported.
    """

    def __init__(self, param_space=1.0):
        self._parse_param_space(param_space)

    @property
    def params(self):
        return np.asarray(self.l)

    @params.setter
    def params(self, theta):
        self.l = theta

    def __call__(self, X, Y=None, eval_gradient=False):
        """ Return the kernel k(X, Y) and optionally its gradient.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Y : array, shape (n_samples_Y, n_features), (optional, default=None)
            Right argument of the returned kernel k(X, Y). If None, k(X, X)
            if evaluated instead.

        eval_gradient : bool (optional, default=False)
            Determines whether the gradient with respect to the kernel
            hyperparameter is determined. Only supported when Y is None.

        Returns
        -------
        K : array, shape (n_samples_X, n_samples_Y)
            Kernel k(X, Y)

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_params)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        if Y is None:
            dists = pdist(X / self.l, metric='sqeuclidean')
            K = np.exp(-.5 * dists)
            # convert from upper-triangular matrix to square matrix
            K = squareform(K)
            np.fill_diagonal(K, 1)
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X / self.l, Y / self.l, metric='sqeuclidean')
            K = np.exp(-.5 * dists)

        if eval_gradient:
            if self.l.shape[0] == 1:
                K_gradient = \
                    (K * squareform(dists) / self.l[0])[:, :, np.newaxis]
                return K, K_gradient
            elif self.l.shape[0] == X.shape[1]:
                # We need to recompute the pairwise dimension-wise distances
                D = (X[:, np.newaxis, :] - X[np.newaxis, :, :]) ** 2 \
                    / (self.l ** 3)
                K_gradient = K[..., np.newaxis] * D
                return K, K_gradient
            else:
                raise Exception("Anisotropic kernels require that the number "
                                "of length scales and features match.")
        else:
            return K


class WhiteKernel(Kernel):
    """ White kernel.

    The main use-case of this kernel is as part of a sum-kernel where it
    explains the noise-component of the signal. Tuning its parameter
    corresponds to estimating the noise-level.
    """

    def __init__(self, param_space=1.0):
        self._parse_param_space(param_space)

    @property
    def params(self):
        return np.asarray([self.c])

    @params.setter
    def params(self, theta):
        self.c = theta[0]

    def __call__(self, X, Y=None, eval_gradient=False):
        """ Return the kernel k(X, Y) and optionally its gradient.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Y : array, shape (n_samples_Y, n_features), (optional, default=None)
            Right argument of the returned kernel k(X, Y). If None, k(X, X)
            if evaluated instead.

        eval_gradient : bool (optional, default=False)
            Determines whether the gradient with respect to the kernel
            hyperparameter is determined. Only supported when Y is None.

        Returns
        -------
        K : array, shape (n_samples_X, n_samples_Y)
            Kernel k(X, Y)

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_params)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        if Y is not None and eval_gradient:
            raise ValueError("Gradient can only be evaluated when Y is None.")

        if Y is None:
            K = self.c * np.eye(X.shape[0])
            if eval_gradient:
                return K, np.eye(X.shape[0])[:, :, np.newaxis]
            else:
                return K
        else:
            K = np.zeros((X.shape[0], Y.shape[0]))
            # entries which are sufficiently similar to be considered identical
            K[cdist(X, Y) < 1e-10] = self.c
            return K


class PairwiseKernel(Kernel):
    """ Wrapper for kernels in sklearn.metrics.pairwise.

    A thin wrapper around the functionality of the kernels in
    sklearn.metrics.pairwise.

    Note: Evaluation of eval_gradient is not analytic but numeric and all
          kernels support only isotropic distances. The parameter gamma is
          specified via the param_space and may be optimized. The other
          kernel parameters are set directly  at initialization and are kept
          fixed.

    Parameters
    ----------
    metric : string, or callable
        The metric to use when calculating kernel between instances in a
        feature array. If metric is a string, it must be one of the metrics
        in pairwise.PAIRWISE_KERNEL_FUNCTIONS.
        If metric is "precomputed", X is assumed to be a kernel matrix.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    `**kwds` : optional keyword parameters
        Any further parameters are passed directly to the kernel function.
    """

    def __init__(self, param_space=1.0, metric="linear", **kwargs):
        self._parse_param_space(param_space)
        self.metric = metric
        self.kwargs = kwargs
        if "gamma" in kwargs:
            raise ValueError(
                "Gamma must not be set directly but via param_space.")

    @property
    def params(self):
        return np.asarray([self.gamma])

    @params.setter
    def params(self, theta):
        self.gamma = theta[0]

    def __call__(self, X, Y=None, eval_gradient=False):
        """ Return the kernel k(X, Y) and optionally its gradient.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Y : array, shape (n_samples_Y, n_features), (optional, default=None)
            Right argument of the returned kernel k(X, Y). If None, k(X, X)
            if evaluated instead.

        eval_gradient : bool (optional, default=False)
            Determines whether the gradient with respect to the kernel
            hyperparameter is determined. Only supported when Y is None.

        Returns
        -------
        K : array, shape (n_samples_X, n_samples_Y)
            Kernel k(X, Y)

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_params)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        K = pairwise_kernels(X, Y, metric=self.metric, gamma=self.gamma,
                             filter_params=True, **self.kwargs)
        if eval_gradient:
            # approximate gradient numerically
            K_gradient = np.empty((K.shape[0], K.shape[1], 1))

            def f(gamma, i, j):  # helper function
                return pairwise_kernels(
                    X, Y, metric=self.metric, gamma=gamma,
                    filter_params=True, **self.kwargs)[i, j]
            # XXX: avoid python for-loops
            for i in range(K.shape[0]):
                for j in range(K.shape[0]):
                    K_gradient[i, j] = \
                        approx_fprime(np.array([self.gamma]),
                                      partial(f, i=i, j=j), 1e-10)

            return K, K_gradient
        else:
            return K
