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
import inspect

import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform

from ..metrics.pairwise import pairwise_kernels
from ..externals import six
from ..base import clone


class Kernel(six.with_metaclass(ABCMeta)):
    """ Base class for all kernels."""

    def get_params(self, deep=True):
        """Get parameters of this kernel.

        Parameters
        ----------
        deep: boolean, optional
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        params = dict()

        # introspect the constructor arguments to find the model parameters
        # to represent
        cls = self.__class__
        init = getattr(cls.__init__, 'deprecated_original', cls.__init__)
        args, varargs, kw, default = inspect.getargspec(init)
        if varargs is not None:
            raise RuntimeError("scikit-learn estimators should always "
                               "specify their parameters in the signature"
                               " of their __init__ (no varargs)."
                               " %s doesn't follow this convention."
                               % (cls, ))
        # Remove 'self' and store remaining arguments in params
        args = args[1:]
        for arg in args:
            params[arg] = getattr(self, arg, None)
        return params

    def clone_with_theta(self, theta):
        """ Returns a clone of self with given hyperparameters theta. """
        cloned = clone(self)
        cloned.theta = theta
        return cloned

    @property
    def n_dims(self):
        """ Returns the number of hyperparameters of the kernel."""
        return self.theta.shape[0]

    @property
    def theta(self):
        theta = []
        for var_name in self.theta_vars:
            if not isinstance(var_name, basestring):  # vector-valued parameter
                var_name, _ = var_name
            theta.append(getattr(self, var_name))
        return np.array(theta)

    @theta.setter
    def theta(self, theta):
        i = 0
        for var_name in self.theta_vars:
            if not isinstance(var_name, basestring):  # vector-valued parameter
                var_name, var_length = var_name
                setattr(self, var_name, theta[i:i + var_length])
                i += var_length
            else:
                setattr(self, var_name, theta[i])
                i += 1

        if i != len(theta):
            raise ValueError("theta has not the correct number of entries."
                             " Should be %d; given are %d"
                             % (i, len(theta)))

    @property
    def bounds(self):
        bounds = []
        for var_name in self.theta_vars:
            if not isinstance(var_name, basestring):  # vector-valued parameter
                var_name, _ = var_name
            bounds.append(getattr(self, var_name + "_bounds"))
        return np.array(bounds)

    @bounds.setter
    def bounds(self, bounds):
        i = 0
        for var_name in self.theta_vars:
            if not isinstance(var_name, basestring):  # vector-valued parameter
                var_name, var_length = var_name
                setattr(self, var_name + "_bounds", bounds[i:i + var_length])
                i += var_length
            else:
                setattr(self, var_name + "_bounds", bounds[i])
                i += 1

        if i != len(bounds):
            raise ValueError("bounds has not the correct number of entries."
                             " Should be %d; given are %d"
                             % (i, len(bounds)))

    def __add__(self, b):
        if not isinstance(b, Kernel):
            return Sum(self, ConstantKernel.from_literal(b))
        return Sum(self, b)

    def __radd__(self, b):
        if not isinstance(b, Kernel):
            return Sum(ConstantKernel.from_literal(b), self)
        return Sum(b, self)

    def __mul__(self, b):
        if not isinstance(b, Kernel):
            return Product(self, ConstantKernel.from_literal(b))
        return Product(self, b)

    def __rmul__(self, b):
        if not isinstance(b, Kernel):
            return Product(ConstantKernel.from_literal(b), self)
        return Product(b, self)

    def __pow__(self, b):
        return Exponentiation(self, b)

    def __eq__(self, b):
        if type(self) != type(b):
            return False
        params_a = self.get_params()
        params_b = b.get_params()
        for key in set(list(params_a.keys()) + list(params_b.keys())):
            if np.any(params_a.get(key, None) != params_b.get(key, None)):
                return False
        return True

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__,
                                 ", ".join(map("{0:.3g}".format, self.theta)))

    @abstractmethod
    def __call__(self, X, Y=None, eval_gradient=False):
        """Evaluate the kernel."""

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        return np.ones(X.shape[0])

    def is_stationary(self):
        """ Returns whether the kernel is stationary. """
        return True


class KernelOperator(Kernel):
    """ Base class for all kernel operators. """

    def __init__(self, k1, k2):
        self.k1 = k1
        self.k2 = k2

    def get_params(self, deep=True):
        """Get parameters of this kernel.

        Parameters
        ----------
        deep: boolean, optional
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        params = dict(k1=self.k1, k2=self.k2)
        return params

    @property
    def theta(self):
        return np.append(self.k1.theta, self.k2.theta)

    @theta.setter
    def theta(self, theta):
        k1_dims = self.k1.n_dims
        self.k1.theta = theta[:k1_dims]
        self.k2.theta = theta[k1_dims:]

    @property
    def bounds(self):
        if self.k1.bounds.size == 0:
            return self.k2.bounds
        if self.k2.bounds.size == 0:
            return self.k1.bounds
        return np.vstack((self.k1.bounds, self.k2.bounds))

    @bounds.setter
    def bounds(self, bounds):
        k1_dims = self.k1.n_dims
        self.k1.bounds = bounds[:k1_dims]
        self.k2.bounds = bounds[k1_dims:]

    def __eq__(self, b):
        return (self.k1 == b.k1 and self.k2 == b.k2) \
            or (self.k1 == b.k2 and self.k2 == b.k1)

    def is_stationary(self):
        """ Returns whether the kernel is stationary. """
        return self.k1.is_stationary() and self.k2.is_stationary()


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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
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

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        return self.k1.diag(X) + self.k2.diag(X)

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
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

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        return self.k1.diag(X) * self.k2.diag(X)

    def __repr__(self):
        return "{0} * {1}".format(self.k1, self.k2)


class Exponentiation(Kernel):
    """ Exponentiate kernel by given exponent. """

    def __init__(self, kernel, exponent):
        self.kernel = kernel
        self.exponent = exponent

    def get_params(self, deep=True):
        """Get parameters of this kernel.

        Parameters
        ----------
        deep: boolean, optional
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        params = dict(kernel=self.kernel, exponent=self.exponent)
        return params

    @property
    def theta(self):
        return self.kernel.theta

    @theta.setter
    def theta(self, theta):
        self.kernel.theta = theta

    @property
    def bounds(self):
        return self.kernel.bounds

    @bounds.setter
    def bounds(self, bounds):
        self.kernel.bounds = bounds

    def __eq__(self, b):
        return (self.kernel == b.kernel and self.exponent == b.exponent)

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        if eval_gradient:
            K, K_gradient = self.kernel(X, Y, eval_gradient=True)
            K_gradient *= \
                self.exponent * K[:, :, np.newaxis] ** (self.exponent - 1)
            return K ** self.exponent, K_gradient
        else:
            K = self.kernel(X, Y, eval_gradient=False)
            return K ** self.exponent

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        return self.kernel.diag(X) ** self.exponent

    def __repr__(self):
        return "{0} ** {1}".format(self.kernel, self.exponent)

    def is_stationary(self):
        """ Returns whether the kernel is stationary. """
        return self.kernel.is_stationary()


class ConstantKernel(Kernel):
    """ Constant kernel.

    Can be used as part of a product-kernel where it scales the magnitude of
    the other factor (kernel) or as part of a sum-kernel, where it modifies
    the mean of the Gaussian process.

    Tunable kernel parameters
    -------------------------
    c : float
        The constant value used for determining the magnitude (product-kernel)
        or offset of mean (sum-kernel).
    """
    def __init__(self, c=1.0, c_bounds=(0, np.inf)):
        self.c = c
        self.c_bounds = c_bounds

        self.theta_vars = ["c"] if c_bounds is not "fixed" else []

    @classmethod
    def from_literal(cls, literal):
        if np.iterable(literal):
            if len(literal) == 1:
                return cls(c=literal[0])
            elif len(literal) == 2:
                return cls(c=(literal[0] + literal[1]) / 2,
                           c_bounds=(literal[0], literal[1]))
            elif len(literal) == 3:
                return cls(c=literal[1], c_bounds=(literal[0], literal[2]))
            else:
                raise ValueError("Cannot interpret literal %s for "
                                 "ConstantKernel." % literal)
        else:
            return cls(literal)

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if Y is None:
            Y = X
        elif eval_gradient:
            raise ValueError("Gradient can only be evaluated when Y is None.")

        K = self.c * np.ones((X.shape[0], Y.shape[0]))
        if eval_gradient:
            if self.c_bounds is not "fixed":
                return K, np.ones((X.shape[0], X.shape[0], 1))
            else:
                return K, np.empty((X.shape[0], X.shape[0], 0))
        else:
            return K

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        return self.c * np.ones(X.shape[0])

    def __repr__(self):
        return "{0:.3g}".format(self.c)


class WhiteKernel(Kernel):
    """ White kernel.

    The main use-case of this kernel is as part of a sum-kernel where it
    explains the noise-component of the signal. Tuning its parameter
    corresponds to estimating the noise-level.

    Tunable kernel parameters
    -------------------------
    c : float
        Parameter controlling the noise level.
    """
    def __init__(self, c=1.0, c_bounds=(0.0, np.inf)):
        self.c = c
        self.c_bounds = c_bounds

        self.theta_vars = ["c"] if c_bounds is not "fixed" else []

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if Y is not None and eval_gradient:
            raise ValueError("Gradient can only be evaluated when Y is None.")

        if Y is None:
            K = self.c * np.eye(X.shape[0])
            if eval_gradient:
                if self.c_bounds is not "fixed":
                    return K, np.eye(X.shape[0])[:, :, np.newaxis]
                else:
                    return K, np.empty((X.shape[0], X.shape[0], 0))
            else:
                return K
        else:
            K = np.zeros((X.shape[0], Y.shape[0]))
            # entries which are sufficiently similar to be considered identical
            K[cdist(X, Y) < 1e-10] = self.c
            return K

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        return self.c * np.ones(X.shape[0])

    def __repr__(self):
        return "{0}(c={1:.3g})".format(self.__class__.__name__, self.c)


class RBF(Kernel):
    """ Radial-basis function kernel (aka squared-exponential kernel).

    Both isotropic and anisotropic version are supported.

    Tunable kernel parameters
    -------------------------
    l : float or array with shape (n_features,), entries > 0
        The length scale of the kernel. If a float, an isotropic kernel is
        used. If an array, an anisotropic kernel is used where each dimension
        of l defines the length-scale of the respective feature dimension.
    """
    def __init__(self, l=1.0, l_bounds=(1e-5, np.inf)):
        if np.iterable(l):
            self.l = np.asarray(l, dtype=np.float)
        else:
            self.l = float(l)
        self.l_bounds = l_bounds

        self.theta_vars = []
        if l_bounds is not "fixed":
            if np.iterable(l):  # anisotropic l needs special care
                self.theta_vars.append(("l", len(l)))
            else:
                self.theta_vars.append("l")

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
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
            if self.l_bounds is "fixed":  # Hyperparameter l kept fixed
                return K, np.empty((X.shape[0], X.shape[0], 0))
            elif not np.iterable(self.l) or self.l.shape[0] == 1:
                K_gradient = \
                    (K * squareform(dists) / self.l)[:, :, np.newaxis]
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

    def __repr__(self):
        if np.iterable(self.l):  # anisotropic
            return "{0}(l=[{1}])".format(self.__class__.__name__,
                                         ", ".join(map("{0:.3g}".format,
                                                   self.l)))
        else:  # isotropic
            return "{0}(l={1:.3g})".format(self.__class__.__name__, self.l)



class RationalQuadratic(Kernel):
    """ Rational Quadratic kernel.

    This kernel can be seen as a scale mixture (an infinite sum) of RBF kernels
    with different characteristic length-scales.

    Only isotropic variant is supported at the moment.

    Tunable kernel parameters
    -------------------------
    alpha : float > 0
        Scale mixture parameter
    l : float > 0
        The length scale of the kernel.
    """
    def __init__(self, alpha=1.0, l=1.0, alpha_bounds=(1e-5, np.inf),
                 l_bounds=(1e-5, np.inf)):
        self.alpha = alpha
        self.l = l
        self.alpha_bounds = alpha_bounds
        self.l_bounds = l_bounds

        self.theta_vars = ["alpha"] if alpha_bounds is not "fixed" else []
        if self.l_bounds is not "fixed":
            self.theta_vars += ["l"]

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if Y is None:
            dists = squareform(pdist(X, metric='sqeuclidean'))
            tmp = dists / (2 * self.alpha * self.l ** 2)
            base = (1 + tmp)
            K = base ** -self.alpha
            np.fill_diagonal(K, 1)
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X, Y, metric='sqeuclidean')
            K = (1 + dists / (2 * self.alpha * self.l ** 2)) ** -self.alpha

        if eval_gradient:
            # gradient with respect to alpha
            if "alpha" in self.theta_vars:
                alpha_gradient = K * (-np.log(base) + tmp / base)
                alpha_gradient = alpha_gradient[:, :, np.newaxis]
            else:  # alpha is kept fixed
                alpha_gradient = np.empty((K.shape[0], K.shape[1], 0))
            # gradient with respect to l
            if "l" in self.theta_vars:
                l_gradient = dists * K / (self.l ** 2 * base)
                l_gradient = l_gradient[:, :, np.newaxis]
            else:  # l is kept fixed
                l_gradient = np.empty((K.shape[0], K.shape[1], 0))

            return K, np.dstack((alpha_gradient, l_gradient))
        else:
            return K

    def __repr__(self):
        return "{0}(alpha={1:.3g}, l={2:.3g})".format(
            self.__class__.__name__, self.alpha, self.l)


class ExpSineSquared(Kernel):
    """ Exp-Sine-Squared kernel.

    This kernel allows modelling periodic functions.

    Only isotropic variant is supported at the moment.

    Tunable kernel parameters
    -------------------------
    l : float > 0
        The length scale of the kernel.
    p : float > 0
        The periodicity of the kernel.
    """

    def __init__(self, l=1.0, p=1.0, l_bounds=(1e-5, np.inf),
                 p_bounds=(1e-5, np.inf)):
        self.l = l
        self.p = p
        self.l_bounds = l_bounds
        self.p_bounds = p_bounds

        self.theta_vars = ["l"] if l_bounds is not "fixed" else []
        if self.p_bounds is not "fixed":
            self.theta_vars += ["p"]

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if Y is None:
            dists = squareform(pdist(X, metric='euclidean'))
            arg = np.pi  * dists / self.p
            sin_of_arg = np.sin(arg)
            K = np.exp(- 2 * (sin_of_arg / self.l) ** 2)
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X, Y, metric='euclidean')
            K = np.exp(- 2 * (np.sin(np.pi / self.p * dists) / self.l) ** 2)

        if eval_gradient:
            cos_of_arg = np.cos(arg)
            # gradient with respect to l
            if "l" in self.theta_vars:
                l_gradient = 4 / self.l**3 * sin_of_arg**2 * K
                l_gradient = l_gradient[:, :, np.newaxis]
            else:  # l is kept fixed
                l_gradient = np.empty((K.shape[0], K.shape[1], 0))
            # gradient with respect to p
            if "p" in self.theta_vars:
                p_gradient = \
                   4 * arg / (self.l**2 * self.p) * cos_of_arg * sin_of_arg * K
                p_gradient = p_gradient[:, :, np.newaxis]
            else:  # p is kept fixed
                p_gradient = np.empty((K.shape[0], K.shape[1], 0))

            return K, np.dstack((l_gradient, p_gradient))
        else:
            return K

    def __repr__(self):
        return "{0}(l={1:.3g}, p={2:.3g})".format(
            self.__class__.__name__, self.l, self.p)


class DotProduct(Kernel):
    """ Dot-Product kernel.

    This kernel is non-stationary.

    Tunable kernel parameters
    -------------------------
    sigma_0 : float >= 0
        Parameter controlling the inhomogenity of the kernel. If sigma_0=0,
        the kernel is homogenous.
    """

    def __init__(self, sigma_0=1.0, sigma_0_bounds=(1e-5, np.inf)):
        self.sigma_0 = sigma_0
        self.sigma_0_bounds = sigma_0_bounds

        self.theta_vars = ["sigma_0"] if sigma_0_bounds is not "fixed" else []

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if Y is None:
            K = np.inner(X, X) + self.sigma_0 ** 2
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            K = np.inner(X, Y) + self.sigma_0 ** 2

        if eval_gradient:
            if self.sigma_0_bounds is not "fixed":
                K_gradient = np.empty((K.shape[0], K.shape[1], 1))
                K_gradient[..., 0] = 2 * self.sigma_0
                return K, K_gradient
            else:
                return K, np.empty((X.shape[0], X.shape[0], 0))
        else:
            return K

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        return (X ** 2).sum(1) + self.sigma_0 ** 2

    def is_stationary(self):
        """ Returns whether the kernel is stationary. """
        return False

    def __repr__(self):
        return "{0}(sigma_0={1:.3g})".format(
            self.__class__.__name__, self.sigma_0)


# adapted from scipy/optimize/optimize.py for functions with 2d output
def _approx_fprime(xk, f, epsilon, args=()):
    f0 = f(*((xk,) + args))
    grad = np.zeros((f0.shape[0], f0.shape[1], len(xk)), float)
    ei = np.zeros((len(xk), ), float)
    for k in range(len(xk)):
        ei[k] = 1.0
        d = epsilon * ei
        grad[:, :, k] = (f(*((xk + d,) + args)) - f0) / d[k]
        ei[k] = 0.0
    return grad


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

    def __init__(self, gamma=1.0, gamma_bounds=(1e-5, np.inf),
                 metric="linear", **kwargs):
        self.gamma = gamma
        self.gamma_bounds = gamma_bounds

        self.theta_vars = ["gamma"] if gamma_bounds is not "fixed" else []

        self.metric = metric
        self.kwargs = kwargs

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

        K_gradient : array (opt.), shape (n_samples_X, n_samples_X, n_dims)
            The gradient of the kernel k(X, X) with repect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        K = pairwise_kernels(X, Y, metric=self.metric, gamma=self.gamma,
                             filter_params=True, **self.kwargs)
        if eval_gradient:
            if self.gamma_bounds is "fixed":
                return K, np.empty((X.shape[0], X.shape[0], 0))
            else:
                # approximate gradient numerically
                def f(gamma):  # helper function
                    return pairwise_kernels(
                        X, Y, metric=self.metric, gamma=gamma,
                        filter_params=True, **self.kwargs)
                return K, _approx_fprime(self.theta, f, 1e-10)
        else:
            return K

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluted more efficiently since only the diagonal is
        evaluated.

        Parameters
        ----------
        X : array, shape (n_samples_X, n_features)
            Left argument of the returned kernel k(X, Y)

        Returns
        -------
        K_diag : array, shape (n_samples_X,)
            Diagonal of kernel k(X, X)
        """
        # We have to fall back to slow way of computing diagonal
        return np.apply_along_axis(self, 1, X)[:, 0]

    def is_stationary(self):
        """ Returns whether the kernel is stationary. """
        return self.metric in ["rbf"]

    def __repr__(self):
        return "{0}(gamma={1}, metric={2})".format(
            self.__class__.__name__, self.gamma, self.metric)
