"""Kernels for Gaussian process regression and classification.

The kernels in this module allow kernel-engineering, i.e., they can be
combined via the "+" and "*" operators or be exponentiated with a scalar
via "**". These sum and product expressions can also contain scalar values,
which are automatically converted to a constant kernel.

All kernels allow (analytic) gradient-based hyperparameter optimization.
The space of hyperparameters can be specified by giving lower und upper
boundaries for the value of each hyperparameter (the search space is thus
rectangular). Instead of specifying bounds, hyperparameters can also be
declared to be "fixed", which causes these hyperparameters to be excluded from
optimization.
"""

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

# Note: this module is strongly inspired by the kernel module of the george
#       package.

from abc import ABCMeta, abstractmethod
from functools import partial
from collections import namedtuple
import inspect
import math

import numpy as np
from scipy.special import kv, gamma
from scipy.spatial.distance import pdist, cdist, squareform

from ..metrics.pairwise import pairwise_kernels
from ..externals import six
from ..base import clone


class Hyperparameter(namedtuple(
     'Hyperparameter',
     ('name', 'value_type', 'bounds', 'n_elements', 'fixed'))):
    # A raw namedtuple is very memory efficient as it packs the attributes
    # in a struct to get rid of the __dict__ of attributes in particular it
    # does not copy the string for the keys on each instance.
    # By deriving a namedtuple class just to introduce the __init__ method we
    # would also reintroduce the __dict__ on the instance. By telling the
    # Python interpreter that this subclass uses static __slots__ instead of
    # dynamic attributes. Furthermore we don't need any additional slot in the
    # subclass so we set __slots__ to the empty tuple.
    __slots__ = ()

    def __new__(cls, name, value_type, bounds, n_elements=1, fixed=None):
        if bounds is not "fixed":
            bounds = np.atleast_2d(bounds)
            if n_elements > 1:  # vector-valued parameter
                if bounds.shape[0] == 1:
                    bounds = np.repeat(bounds, n_elements, 0)
                elif bounds.shape[0] != n_elements:
                    raise ValueError("Bounds on %s should have either 1 or "
                                     "%d dimensions. Given are %d"
                                     % (name, n_elements, bounds.shape[0]))

        if fixed is None:
             fixed = bounds is "fixed"
        return super(Hyperparameter, cls).__new__(
            cls, name, value_type, bounds, n_elements, fixed)


class Kernel(six.with_metaclass(ABCMeta)):
    """Base class for all kernels."""

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
            raise RuntimeError("scikit-learn kernels should always "
                               "specify their parameters in the signature"
                               " of their __init__ (no varargs)."
                               " %s doesn't follow this convention."
                               % (cls, ))
        # Remove 'self' and store remaining arguments in params
        args = args[1:]
        for arg in args:
            params[arg] = getattr(self, arg, None)
        return params

    def set_params(self, **params):
        """Set the parameters of this kernel.

        The method works on simple kernels as well as on nested kernels.
        The latter have parameters of the form ``<component>__<parameter>``
        so that it's possible to update each component of a nested object.

        Returns
        -------
        self
        """
        if not params:
            # Simple optimisation to gain speed (inspect is slow)
            return self
        valid_params = self.get_params(deep=True)
        for key, value in six.iteritems(params):
            split = key.split('__', 1)
            if len(split) > 1:
                # nested objects case
                name, sub_name = split
                if name not in valid_params:
                    raise ValueError('Invalid parameter %s for kernel %s. '
                                     'Check the list of available parameters '
                                     'with `kernel.get_params().keys()`.' %
                                     (name, self))
                sub_object = valid_params[name]
                sub_object.set_params(**{sub_name: value})
            else:
                # simple objects case
                if key not in valid_params:
                    raise ValueError('Invalid parameter %s for kernel %s. '
                                     'Check the list of available parameters '
                                     'with `kernel.get_params().keys()`.' %
                                     (key, self.__class__.__name__))
                setattr(self, key, value)
        return self

    def clone_with_theta(self, theta):
        """Returns a clone of self with given hyperparameters theta. """
        cloned = clone(self)
        cloned.theta = theta
        return cloned

    @property
    def n_dims(self):
        """Returns the number of non-fixed hyperparameters of the kernel."""
        return self.theta.shape[0]

    @property
    def hyperparameters(self):
        """Returns a list of all hyperparameter."""
        r = []
        for attr, value in self.__dict__.items():
            if attr.startswith("hyperparameter_"):
                r.append(value)
        return r

    @property
    def theta(self):
        """Returns the (flattened, log-transformed) non-fixed hyperparameters.

        Note that theta are typically the log-transformed values of the
        kernel's hyperparameters as this representation of the search space
        is more amenable for hyperparameter search, as hyperparameters like
        length-scales naturally live on a log-scale.

        Returns
        -------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        theta = []
        for hyperparameter in self.hyperparameters:
            if not hyperparameter.fixed:
                theta.append(getattr(self, hyperparameter.name))
        if len(theta) > 0:
            return np.log(np.hstack(theta))
        else:
            return np.array([])

    @theta.setter
    def theta(self, theta):
        """Sets the (flattened, log-transformed) non-fixed hyperparameters.

        Parameters
        ----------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        i = 0
        for hyperparameter in self.hyperparameters:
            if hyperparameter.fixed:
                continue
            if hyperparameter.n_elements > 1:
                # vector-valued parameter
                setattr(self, hyperparameter.name,
                        np.exp(theta[i:i + hyperparameter.n_elements]))
                i += hyperparameter.n_elements
            else:
                setattr(self, hyperparameter.name, np.exp(theta[i]))
                i += 1

        if i != len(theta):
            raise ValueError("theta has not the correct number of entries."
                             " Should be %d; given are %d"
                             % (i, len(theta)))

    @property
    def bounds(self):
        """Returns the bounds on the kernel's hyperparameters theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        bounds = []
        for hyperparameter in self.hyperparameters:
            if not hyperparameter.fixed:
                bounds.append(hyperparameter.bounds)
        if len(bounds) > 0:
            return np.log(np.vstack(bounds))
        else:
            return np.array([])

    @bounds.setter
    def bounds(self, bounds):
        """Sets the bounds on the kernel's hyperparameters theta.

        Parameters
        ----------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        bounds_exp = np.exp(bounds)
        i = 0
        for hyperparameter in self.hyperparameters:
            if hyperparameter.n_elements > 1:  # vector-valued parameter
                setattr(self, "hyperparameter_" + hyperparameter.name,
                        Hyperparameter(
                            hyperparameter.name, hyperparameter.value_type,
                            bounds_exp[i:i + hyperparameter.n_elements],
                            hyperparameter.n_elements))
                setattr(self, hyperparameter.name + "_bounds",
                        bounds_exp[i:i + hyperparameter.n_elements])
                i += var_length
            else:
                setattr(self, "hyperparameter_" + hyperparameter.name,
                        Hyperparameter(hyperparameter.name,
                                       hyperparameter.value_type,
                                       bounds_exp[i]))
                setattr(self, hyperparameter.name + "_bounds", bounds_exp[i])
                i += 1

        if i != len(bounds):
            raise ValueError("bounds has not the correct number of entries."
                             " Should be %d; given are %d"
                             % (i, len(bounds)))

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
        it can be evaluated more efficiently since only the diagonal is
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
        """Returns whether the kernel is stationary. """
        return True


class CompoundKernel(Kernel):
    """Kernel which is composed of a set of other kernels."""

    def __init__(self, kernels):
        self.kernels = kernels

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
        return dict(kernels=kernels)

    @property
    def theta(self):
        """Returns the (flattened, log-transformed) non-fixed hyperparameters.

        Note that theta are typically the log-transformed values of the
        kernel's hyperparameters as this representation of the search space
        is more amenable for hyperparameter search, as hyperparameters like
        length-scales naturally live on a log-scale.

        Returns
        -------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        return np.hstack([kernel.theta for kernel in self.kernels])

    @theta.setter
    def theta(self, theta):
        """Sets the (flattened, log-transformed) non-fixed hyperparameters.

        Parameters
        ----------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        k_dims = self.k1.n_dims
        for i, kernel in enumerate(self.kernels):
            kernel.theta = theta[i*k_dims:(i+1)*k_dims]

    @property
    def bounds(self):
        """Returns the bounds on the kernel's hyperparameters theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        return np.vstack([kernel.bounds for kernel in self.kernels])

    @bounds.setter
    def bounds(self, bounds):
        """Sets the bounds on the kernel's hyperparameters theta.

        Parameters
        ----------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        k1_dims = self.k1.n_dims
        for i, kernel in enumerate(self.kernels):
            kernel.bounds = bounds[i*k_dims:(i+1)*k_dims]

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

        Note that this compound kernel returns the results of all simple kernel
        stacked along an additional axis.

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
        K : array, shape (n_samples_X, n_samples_Y, n_kernels)
            Kernel k(X, Y)

        K_gradient : array, shape (n_samples_X, n_samples_X, n_dims, n_kernels)
            The gradient of the kernel k(X, X) with respect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        if eval_gradient:
            K = []
            K_grad = []
            for kernel in self.kernels:
                K_single, K_grad_single = kernel(X, Y, eval_gradient)
                K.append(K_single)
                K_grad.append(K_grad_single[..., np.newaxis])
            return np.dstack(K), np.concatenate(K_grad, 3)
        else:
            return np.dstack([kernel(X, Y, eval_gradient)
                              for kernel in self.kernels])

    def __eq__(self, b):
        if type(self) != type(b) or len(self.kernels) != len(b.kernels):
            return False
        return np.all([self.kernels[i] == b.kernels[i]
                       for i in range(len(self.kernels))])

    def is_stationary(self):
        """Returns whether the kernel is stationary. """
        return np.all([kernel.is_stationary() for kernel in self.kernels])


class KernelOperator(Kernel):
    """Base class for all kernel operators. """

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
        if deep:
            deep_items = self.k1.get_params().items()
            params.update(('k1__' + k, val) for k, val in deep_items)
            deep_items = self.k2.get_params().items()
            params.update(('k2__' + k, val) for k, val in deep_items)

        return params

    @property
    def hyperparameters(self):
        """Returns a list of all hyperparameter."""
        r = []
        for hyperparameter in self.k1.hyperparameters:
            r.append(Hyperparameter("k1__" + hyperparameter.name,
                                    hyperparameter.value_type,
                                    hyperparameter.bounds,
                                    hyperparameter.n_elements))
        for hyperparameter in self.k2.hyperparameters:
            r.append(Hyperparameter("k2__" + hyperparameter.name,
                                    hyperparameter.value_type,
                                    hyperparameter.bounds,
                                    hyperparameter.n_elements))
        return r

    @property
    def theta(self):
        """Returns the (flattened, log-transformed) non-fixed hyperparameters.

        Note that theta are typically the log-transformed values of the
        kernel's hyperparameters as this representation of the search space
        is more amenable for hyperparameter search, as hyperparameters like
        length-scales naturally live on a log-scale.

        Returns
        -------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        return np.append(self.k1.theta, self.k2.theta)

    @theta.setter
    def theta(self, theta):
        """Sets the (flattened, log-transformed) non-fixed hyperparameters.

        Parameters
        ----------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        k1_dims = self.k1.n_dims
        self.k1.theta = theta[:k1_dims]
        self.k2.theta = theta[k1_dims:]

    @property
    def bounds(self):
        """Returns the bounds on the kernel's hyperparameters theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        if self.k1.bounds.size == 0:
            return self.k2.bounds
        if self.k2.bounds.size == 0:
            return self.k1.bounds
        return np.vstack((self.k1.bounds, self.k2.bounds))

    @bounds.setter
    def bounds(self, bounds):
        """Sets the bounds on the kernel's hyperparameters theta.

        Parameters
        ----------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        k1_dims = self.k1.n_dims
        self.k1.bounds = bounds[:k1_dims]
        self.k2.bounds = bounds[k1_dims:]

    def __eq__(self, b):
        if type(self) != type(b):
            return False
        return (self.k1 == b.k1 and self.k2 == b.k2) \
            or (self.k1 == b.k2 and self.k2 == b.k1)

    def is_stationary(self):
        """Returns whether the kernel is stationary. """
        return self.k1.is_stationary() and self.k2.is_stationary()


class Sum(KernelOperator):
    """Sum-kernel k1 + k2 of two kernels k1 and k2.

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
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
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
        it can be evaluated more efficiently since only the diagonal is
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
    """Product-kernel k1 * k2 of two kernels k1 and k2.

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
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
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
        it can be evaluated more efficiently since only the diagonal is
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
    """Exponentiate kernel by given exponent.

    The resulting kernel is defined as
    k_exp(X, Y) = k(X, Y) ** exponent

    Parameters
    ----------
    kernel : Kernel object
        The base kernel

    exponent : float
        The exponent for the base kernel

    """
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
        if deep:
            deep_items = self.kernel.get_params().items()
            params.update(('kernel__' + k, val) for k, val in deep_items)
        return params

    @property
    def hyperparameters(self):
        """Returns a list of all hyperparameter."""
        r = []
        for hyperparameter in self.kernel.hyperparameters:
            r.append(Hyperparameter("kernel__" + hyperparameter.name,
                                    hyperparameter.value_type,
                                    hyperparameter.bounds,
                                    hyperparameter.n_elements))
        return r

    @property
    def theta(self):
        """Returns the (flattened, log-transformed) non-fixed hyperparameters.

        Note that theta are typically the log-transformed values of the
        kernel's hyperparameters as this representation of the search space
        is more amenable for hyperparameter search, as hyperparameters like
        length-scales naturally live on a log-scale.

        Returns
        -------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        return self.kernel.theta

    @theta.setter
    def theta(self, theta):
        """Sets the (flattened, log-transformed) non-fixed hyperparameters.

        Parameters
        ----------
        theta : array, shape (n_dims,)
            The non-fixed, log-transformed hyperparameters of the kernel
        """
        self.kernel.theta = theta

    @property
    def bounds(self):
        """Returns the bounds on the kernel's hyperparameters theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        return self.kernel.bounds

    @bounds.setter
    def bounds(self, bounds):
        """Sets the bounds on the kernel's hyperparameters theta.

        Parameters
        ----------
        bounds : array, shape (n_dims, 2)
            The bounds on the kernel's hyperparameters theta
        """
        self.kernel.bounds = bounds

    def __eq__(self, b):
        if type(self) != type(b):
            return False
        return (self.kernel == b.kernel and self.exponent == b.exponent)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
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
        it can be evaluated more efficiently since only the diagonal is
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
        """Returns whether the kernel is stationary. """
        return self.kernel.is_stationary()


class ConstantKernel(Kernel):
    """Constant kernel.

    Can be used as part of a product-kernel where it scales the magnitude of
    the other factor (kernel) or as part of a sum-kernel, where it modifies
    the mean of the Gaussian process.

    k(x_1, x_2) = c for all x_1, x_2

    Parameters
    ----------
    c : float, default: 1.0
        The constant value which defines the covariance: k(x_1, x_2) = c

    c_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on c
    """
    def __init__(self, c=1.0, c_bounds=(1e-5, 1e5)):
        self.c = c
        self.c_bounds = c_bounds

        self.hyperparameter_c = Hyperparameter("c", "numeric", c_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
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
            if not self.hyperparameter_c.fixed:
                return K, self.c * np.ones((X.shape[0], X.shape[0], 1))
            else:
                return K, np.empty((X.shape[0], X.shape[0], 0))
        else:
            return K

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluated more efficiently since only the diagonal is
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
        return "{0:.3g}**2".format(np.sqrt(self.c))


class WhiteKernel(Kernel):
    """White kernel.

    The main use-case of this kernel is as part of a sum-kernel where it
    explains the noise-component of the signal. Tuning its parameter
    corresponds to estimating the noise-level.

    k(x_1, x_2) = c if x_1 == x_2 else 0

    Parameters
    ----------
    c : float, default: 1.0
        Parameter controlling the noise level

    c_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on c
    """
    def __init__(self, c=1.0, c_bounds=(1e-5, 1e5)):
        self.c = c
        self.c_bounds = c_bounds

        self.hyperparameter_c = Hyperparameter("c", "numeric", c_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if Y is not None and eval_gradient:
            raise ValueError("Gradient can only be evaluated when Y is None.")

        if Y is None:
            K = self.c * np.eye(X.shape[0])
            if eval_gradient:
                if not self.hyperparameter_c.fixed:
                    return K, self.c * np.eye(X.shape[0])[:, :, np.newaxis]
                else:
                    return K, np.empty((X.shape[0], X.shape[0], 0))
            else:
                return K
        else:
            return np.zeros((X.shape[0], Y.shape[0]))

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluated more efficiently since only the diagonal is
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
    """Radial-basis function kernel (aka squared-exponential kernel).

    The RBF kernel is a stationary kernel. It is also known as the
    "squared exponential" kernel. It is parameterized by a length-scale
    parameter l>0, which can either be a scalar (isotropic variant of
    the kernel) or a vector with the same number of dimensions as the inputs
    X (anisotropic variant of the kernel). The kernel given by:

    k(x_i, x_j) = exp(-1 / 2 d(x_i / l, x_j / l)^2)

    This kernel is infinitely differentiable, which implies that GPs with this
    kernel as covariance function have mean square derivatives of all orders,
    and are thus very smooth.

    Parameters
    -----------
    l : float or array with shape (n_features,), entries > 0, default: 1.0
        The length scale of the kernel. If a float, an isotropic kernel is
        used. If an array, an anisotropic kernel is used where each dimension
        of l defines the length-scale of the respective feature dimension.

    l_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on l
    """
    def __init__(self, l=1.0, l_bounds=(1e-5, 1e5)):
        if np.iterable(l):
            if len(l) > 1:
                self.anisotropic = True
                self.l = np.asarray(l, dtype=np.float)
            else:
                self.anisotropic = False
                self.l = float(l[0])
        else:
            self.anisotropic = False
            self.l = float(l)
        self.l_bounds = l_bounds

        if self.anisotropic:  # anisotropic l
            self.hyperparameter_l = \
                Hyperparameter("l", "numeric", l_bounds, len(l))
        else:
            self.hyperparameter_l = \
                Hyperparameter("l", "numeric", l_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if self.anisotropic and X.shape[1] != self.l.shape[0]:
            raise Exception("Anisotropic kernel must have the same number of "
                            "dimensions as data (%d!=%d)"
                            % (self.l.shape[0], X.shape[1]))

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
            if self.hyperparameter_l.fixed:
                # Hyperparameter l kept fixed
                return K, np.empty((X.shape[0], X.shape[0], 0))
            elif not self.anisotropic or self.l.shape[0] == 1:
                K_gradient = \
                    (K * squareform(dists))[:, :, np.newaxis]
                return K, K_gradient
            elif self.anisotropic:
                # We need to recompute the pairwise dimension-wise distances
                K_gradient = (X[:, np.newaxis, :] - X[np.newaxis, :, :]) ** 2 \
                    / (self.l ** 2)
                K_gradient *= K[..., np.newaxis]
                return K, K_gradient
            else:
                raise Exception("Anisotropic kernels require that the number "
                                "of length scales and features match.")
        else:
            return K

    def __repr__(self):
        if self.anisotropic:
            return "{0}(l=[{1}])".format(self.__class__.__name__,
                                         ", ".join(map("{0:.3g}".format,
                                                   self.l)))
        else:  # isotropic
            return "{0}(l={1:.3g})".format(self.__class__.__name__, self.l)


class Matern(RBF):
    """ Matern kernel.

    The class of Matern kernels is a generalization of the RBF and the
    absolute exponential kernel parameterized by an additional parameter
    nu. The smaller nu, the less smooth the approximated function is.
    For nu=inf, the kernel becomes equivalent to the RBF kernel and for nu=0.5
    to the absolute exponential kernel. Important intermediate values are
    nu=1.5 (once differentiable functions) and nu=2.5 (twice differentiable
    functions).

    See Rasmussen and Williams 2006, pp84 for details regarding the
    different variants of the Matern kernel.

    Parameters
    -----------
    l : float or array with shape (n_features,), entries > 0, default: 1.0
        The length scale of the kernel. If a float, an isotropic kernel is
        used. If an array, an anisotropic kernel is used where each dimension
        of l defines the length-scale of the respective feature dimension.

    l_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on l

    nu: float, default: 1.5
        The parameter nu controlling the smoothness of the learned function.
        The smaller nu, the less smooth the approximated function is.
        For nu=inf, the kernel becomes equivalent to the RBF kernel and for
        nu=0.5 to the absolute exponential kernel. Important intermediate
        values are nu=1.5 (once differentiable functions) and nu=2.5
        (twice differentiable functions). Note that values of nu not in
        [0.5, 1.5, 2.5, inf] incur a considerably higher computational cost
        (appr. 10 times higher) since they require to evaluate the modified
        Bessel function. Furthermore, in contrast to l, nu is kept fixed to
        its initial value and not optimized.
    """
    def __init__(self, l=1.0, l_bounds=(1e-5, 1e5), nu=1.5):
        super(Matern, self).__init__(l, l_bounds)
        self.nu = nu

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if self.anisotropic and X.shape[1] != self.l.shape[0]:
            raise Exception("Anisotropic kernel must have the same number of "
                            "dimensions as data (%d!=%d)"
                            % (self.l.shape[0], X.shape[1]))

        if Y is None:
            dists = pdist(X / self.l, metric='euclidean')
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X / self.l, Y / self.l, metric='euclidean')

        if self.nu == 0.5:
            K = np.exp(-dists)
        elif self.nu == 1.5:
            K = dists * math.sqrt(3)
            K = (1. + K) * np.exp(-K)
        elif self.nu == 2.5:
            K = dists * math.sqrt(5)
            K = (1. + K + K ** 2 / 3.0) * np.exp(-K)
        else:  # general case; expensive to evaluate
            K = dists
            K[K == 0.0] += np.finfo(float).eps  # strict zeros result in nan
            tmp = (math.sqrt(2 * self.nu) * K)
            K.fill((2 ** (1. - self.nu)) / gamma(self.nu))
            K *= tmp ** self.nu
            K *= kv(self.nu, tmp)

        if Y is None:
            # convert from upper-triangular matrix to square matrix
            K = squareform(K)
            np.fill_diagonal(K, 1)

        if eval_gradient:
            if self.hyperparameter_l.fixed:
                # Hyperparameter l kept fixed
                K_gradient = np.empty((X.shape[0], X.shape[0], 0))
                return K, K_gradient

            # We need to recompute the pairwise dimension-wise distances
            if self.anisotropic:
                D = (X[:, np.newaxis, :] - X[np.newaxis, :, :])**2 \
                    / (self.l ** 2)
            else:
                D = squareform(dists**2)[:, :, np.newaxis]

            if self.nu == 0.5:
                K_gradient = K[..., np.newaxis] * D \
                    / np.sqrt(D.sum(2))[:, :, np.newaxis]
                K_gradient[~np.isfinite(K_gradient)] = 0
            elif self.nu == 1.5:
                K_gradient = \
                    3 * D * np.exp(-np.sqrt(3 * D.sum(-1)))[..., np.newaxis]
            elif self.nu == 2.5:
                tmp = np.sqrt(5 * D.sum(-1))[..., np.newaxis]
                K_gradient = 5.0/3.0 * D * (tmp + 1) * np.exp(-tmp)
            else:
                # approximate gradient numerically
                def f(theta):  # helper function
                    return self.clone_with_theta(theta)(X, Y)
                return K, _approx_fprime(self.theta, f, 1e-10)

            if not self.anisotropic:
                return K, K_gradient[:, :].sum(-1)[:, :, np.newaxis]
            else:
                return K, K_gradient
        else:
            return K

    def __repr__(self):
        if self.anisotropic:
            return "{0}(l=[{1}], nu={2:.3g})".format(
                self.__class__.__name__,
                ", ".join(map("{0:.3g}".format, self.l)),
                self.nu)
        else:  # isotropic
            return "{0}(l={1:.3g}, nu={2:.3g})".format(
                self.__class__.__name__, self.l, self.nu)


class RationalQuadratic(Kernel):
    """Rational Quadratic kernel.

    The RationalQuadratic kernel can be seen as a scale mixture (an infinite
    sum) of RBF kernels with different characteristic length-scales. It is
    parameterized by a length-scale parameter l>0 and a scale mixture parameter
    alpha>0 Only the isotropic variant where l is a scalar is supported at the
    moment. The kernel given by:

    k(x_i, x_j) = (1 + d(x_i, x_j)^2 / (2*alpha l^2))^-alpha

    Parameters
    ----------
    l : float > 0, default: 1.0
        The length scale of the kernel.

    alpha : float > 0, default: 1.0
        Scale mixture parameter

    l_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on l

    alpha_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on alpha
    """
    def __init__(self, l=1.0, alpha=1.0, l_bounds=(1e-5, 1e5),
                 alpha_bounds=(1e-5, 1e5)):
        self.l = l
        self.alpha = alpha
        self.l_bounds = l_bounds
        self.alpha_bounds = alpha_bounds

        self.hyperparameter_l = Hyperparameter("l", "numeric", l_bounds)
        self.hyperparameter_alpha = \
             Hyperparameter("alpha", "numeric", alpha_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
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
            # gradient with respect to l
            if not self.hyperparameter_l.fixed:
                l_gradient = dists * K / (self.l ** 2 * base)
                l_gradient = l_gradient[:, :, np.newaxis]
            else:  # l is kept fixed
                l_gradient = np.empty((K.shape[0], K.shape[1], 0))

            # gradient with respect to alpha
            if not self.hyperparameter_alpha.fixed:
                alpha_gradient = \
                    K * (-self.alpha * np.log(base)
                         + dists / (2 * self.l ** 2 * base))
                alpha_gradient = alpha_gradient[:, :, np.newaxis]
            else:  # alpha is kept fixed
                alpha_gradient = np.empty((K.shape[0], K.shape[1], 0))

            return K, np.dstack((l_gradient, alpha_gradient))
        else:
            return K

    def __repr__(self):
        return "{0}(alpha={1:.3g}, l={2:.3g})".format(
            self.__class__.__name__, self.alpha, self.l)


class ExpSineSquared(Kernel):
    """Exp-Sine-Squared kernel.

    The ExpSineSquared kernel allows modeling periodic functions. It is
    parameterized by a length-scale parameter l>0 and a periodicity parameter
    p>0. Only the isotropic variant where l is a scalar is supported at the
    moment. The kernel given by:

    k(x_i, x_j) =  exp(-2 sin(\pi / p * d(x_i, x_j)) / l)^2

    Parameters
    ----------
    l : float > 0, default: 1.0
        The length scale of the kernel.

    p : float > 0, default: 1.0
        The periodicity of the kernel.

    l_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on l

    p_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on p
    """
    def __init__(self, l=1.0, p=1.0, l_bounds=(1e-5, 1e5),
                 p_bounds=(1e-5, 1e5)):
        self.l = l
        self.p = p
        self.l_bounds = l_bounds
        self.p_bounds = p_bounds

        self.hyperparameter_l = Hyperparameter("l", "numeric", l_bounds)
        self.hyperparameter_p = Hyperparameter("p", "numeric", p_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        if Y is None:
            dists = squareform(pdist(X, metric='euclidean'))
            arg = np.pi * dists / self.p
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
            if not self.hyperparameter_l.fixed:
                l_gradient = 4 / self.l**2 * sin_of_arg**2 * K
                l_gradient = l_gradient[:, :, np.newaxis]
            else:  # l is kept fixed
                l_gradient = np.empty((K.shape[0], K.shape[1], 0))
            # gradient with respect to p
            if not self.hyperparameter_p.fixed:
                p_gradient = 4 * arg / self.l**2 * cos_of_arg \
                    * sin_of_arg * K
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
    """Dot-Product kernel.

    The DotProduct kernel is non-stationary and can be obtained from linear
    regression by putting N(0, 1) priors on the coefficients of x_d (d = 1, . .
    . , D) and a prior of N(0, \sigma_0^2) on the bias. The DotProduct kernel
    is invariant to a rotation of the coordinates about the origin, but not
    translations. It is parameterized by a parameter sigma_0^2. For
    sigma_0^2 =0, the kernel is called the homogeneous linear kernel, otherwise
    it is inhomogeneous. The kernel is given by

    k(x_i, x_j) = sigma_0 ^ 2 + x_i \cdot x_j

    The DotProduct kernel is commonly combined with exponentiation.

    Parameters
    ----------
    sigma_0 : float >= 0, default: 1.0
        Parameter controlling the inhomogenity of the kernel. If sigma_0=0,
        the kernel is homogenous.

    sigma_0_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on l
    """

    def __init__(self, sigma_0=1.0, sigma_0_bounds=(1e-5, 1e5)):
        self.sigma_0 = sigma_0
        self.sigma_0_bounds = sigma_0_bounds

        self.hyperparameter_sigma_0 = \
            Hyperparameter("sigma_0", "numeric", sigma_0_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
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
            if not self.hyperparameter_sigma_0.fixed:
                K_gradient = np.empty((K.shape[0], K.shape[1], 1))
                K_gradient[..., 0] = 2 * self.sigma_0 ** 2
                return K, K_gradient
            else:
                return K, np.empty((X.shape[0], X.shape[0], 0))
        else:
            return K

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluated more efficiently since only the diagonal is
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
        return np.einsum('ij,ij->i', X, X) + self.sigma_0 ** 2

    def is_stationary(self):
        """Returns whether the kernel is stationary. """
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
    """Wrapper for kernels in sklearn.metrics.pairwise.

    A thin wrapper around the functionality of the kernels in
    sklearn.metrics.pairwise.

    Note: Evaluation of eval_gradient is not analytic but numeric and all
          kernels support only isotropic distances. The parameter gamma is
          specified via the param_space and may be optimized. The other
          kernel parameters are set directly  at initialization and are kept
          fixed.

    Parameters
    ----------
    gamma: float >= 0, default: 1.0
        Parameter gamma of the pairwise kernel specified by metric

    gamma_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on gamma

    metric : string, or callable, default: "linear"
        The metric to use when calculating kernel between instances in a
        feature array. If metric is a string, it must be one of the metrics
        in pairwise.PAIRWISE_KERNEL_FUNCTIONS.
        If metric is "precomputed", X is assumed to be a kernel matrix.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    pairwise_kernels_kwargs : dict, default: None
        All entries of this dict (if any) are passed as keyword arguments to
        the pairwise kernel function.
    """

    def __init__(self, gamma=1.0, gamma_bounds=(1e-5, 1e5), metric="linear",
                 pairwise_kernels_kwargs=None):
        self.gamma = gamma
        self.gamma_bounds = gamma_bounds

        self.hyperparameter_gamma = \
            Hyperparameter("gamma", "numeric", gamma_bounds)

        self.metric = metric
        if pairwise_kernels_kwargs is not None:
            self.pairwise_kernels_kwargs = pairwise_kernels_kwargs
        else:
            self.pairwise_kernels_kwargs = {}

    def __call__(self, X, Y=None, eval_gradient=False):
        """Return the kernel k(X, Y) and optionally its gradient.

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
            The gradient of the kernel k(X, X) with respect to the
            hyperparameter of the kernel. Only returned when eval_gradient
            is True.
        """
        X = np.atleast_2d(X)
        K = pairwise_kernels(X, Y, metric=self.metric, gamma=self.gamma,
                             filter_params=True,
                             **self.pairwise_kernels_kwargs)
        if eval_gradient:
            if self.hyperparameter_gamma.fixed:
                return K, np.empty((X.shape[0], X.shape[0], 0))
            else:
                # approximate gradient numerically
                def f(gamma):  # helper function
                    return pairwise_kernels(
                        X, Y, metric=self.metric, gamma=np.exp(gamma),
                        filter_params=True, **self.pairwise_kernels_kwargs)
                return K, _approx_fprime(self.theta, f, 1e-10)
        else:
            return K

    def diag(self, X):
        """Returns the diagonal of the kernel k(X, X).

        The result of this method is identical to np.diag(self(X)); however,
        it can be evaluated more efficiently since only the diagonal is
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
        """Returns whether the kernel is stationary. """
        return self.metric in ["rbf"]

    def __repr__(self):
        return "{0}(gamma={1}, metric={2})".format(
            self.__class__.__name__, self.gamma, self.metric)
