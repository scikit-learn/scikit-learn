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
# License: BSD 3 clause

# Note: this module is strongly inspired by the kernel module of the george
#       package.

from abc import ABCMeta, abstractmethod
from collections import namedtuple
import math

import numpy as np
from scipy.special import kv, gamma
from scipy.spatial.distance import pdist, cdist, squareform

from ..metrics.pairwise import pairwise_kernels
from ..externals import six
from ..base import clone
from sklearn.externals.funcsigs import signature


def _check_length_scale(X, length_scale):
    length_scale = np.squeeze(length_scale).astype(float)
    if np.ndim(length_scale) > 1:
        raise ValueError("length_scale cannot be of dimension greater than 1")
    if np.ndim(length_scale) == 1 and X.shape[1] != length_scale.shape[0]:
        raise ValueError("Anisotropic kernel must have the same number of "
                         "dimensions as data (%d!=%d)"
                         % (length_scale.shape[0], X.shape[1]))
    return length_scale


class Hyperparameter(namedtuple('Hyperparameter',
                                ('name', 'value_type', 'bounds',
                                 'n_elements', 'fixed'))):
    """A kernel hyperparameter's specification in form of a namedtuple.

    .. versionadded:: 0.18

    Attributes
    ----------
    name : string
        The name of the hyperparameter. Note that a kernel using a
        hyperparameter with name "x" must have the attributes self.x and
        self.x_bounds

    value_type : string
        The type of the hyperparameter. Currently, only "numeric"
        hyperparameters are supported.

    bounds : pair of floats >= 0 or "fixed"
        The lower and upper bound on the parameter. If n_elements>1, a pair
        of 1d array with n_elements each may be given alternatively. If
        the string "fixed" is passed as bounds, the hyperparameter's value
        cannot be changed.

    n_elements : int, default=1
        The number of elements of the hyperparameter value. Defaults to 1,
        which corresponds to a scalar hyperparameter. n_elements > 1
        corresponds to a hyperparameter which is vector-valued,
        such as, e.g., anisotropic length-scales.

    fixed : bool, default: None
        Whether the value of this hyperparameter is fixed, i.e., cannot be
        changed during hyperparameter tuning. If None is passed, the "fixed" is
        derived based on the given bounds.

    """
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
        if not isinstance(bounds, six.string_types) or bounds != "fixed":
            bounds = np.atleast_2d(bounds)
            if n_elements > 1:  # vector-valued parameter
                if bounds.shape[0] == 1:
                    bounds = np.repeat(bounds, n_elements, 0)
                elif bounds.shape[0] != n_elements:
                    raise ValueError("Bounds on %s should have either 1 or "
                                     "%d dimensions. Given are %d"
                                     % (name, n_elements, bounds.shape[0]))

        if fixed is None:
            fixed = isinstance(bounds, six.string_types) and bounds == "fixed"
        return super(Hyperparameter, cls).__new__(
            cls, name, value_type, bounds, n_elements, fixed)

    # This is mainly a testing utility to check that two hyperparameters
    # are equal.
    def __eq__(self, other):
        return (self.name == other.name and
                self.value_type == other.value_type and
                np.all(self.bounds == other.bounds) and
                self.n_elements == other.n_elements and
                self.fixed == other.fixed)


class Kernel(six.with_metaclass(ABCMeta)):
    """Base class for all kernels.

    .. versionadded:: 0.18
    """

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
        init_sign = signature(init)
        args, varargs = [], []
        for parameter in init_sign.parameters.values():
            if (parameter.kind != parameter.VAR_KEYWORD and
                    parameter.name != 'self'):
                args.append(parameter.name)
            if parameter.kind == parameter.VAR_POSITIONAL:
                varargs.append(parameter.name)

        if len(varargs) != 0:
            raise RuntimeError("scikit-learn kernels should always "
                               "specify their parameters in the signature"
                               " of their __init__ (no varargs)."
                               " %s doesn't follow this convention."
                               % (cls, ))
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
        """Returns a list of all hyperparameter specifications."""
        r = []
        for attr in dir(self):
            if attr.startswith("hyperparameter_"):
                r.append(getattr(self, attr))
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
        params = self.get_params()
        for hyperparameter in self.hyperparameters:
            if not hyperparameter.fixed:
                theta.append(params[hyperparameter.name])
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
        params = self.get_params()
        i = 0
        for hyperparameter in self.hyperparameters:
            if hyperparameter.fixed:
                continue
            if hyperparameter.n_elements > 1:
                # vector-valued parameter
                params[hyperparameter.name] = np.exp(
                    theta[i:i + hyperparameter.n_elements])
                i += hyperparameter.n_elements
            else:
                params[hyperparameter.name] = np.exp(theta[i])
                i += 1

        if i != len(theta):
            raise ValueError("theta has not the correct number of entries."
                             " Should be %d; given are %d"
                             % (i, len(theta)))
        self.set_params(**params)

    @property
    def bounds(self):
        """Returns the log-transformed bounds on the theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The log-transformed bounds on the kernel's hyperparameters theta
        """
        bounds = []
        for hyperparameter in self.hyperparameters:
            if not hyperparameter.fixed:
                bounds.append(hyperparameter.bounds)
        if len(bounds) > 0:
            return np.log(np.vstack(bounds))
        else:
            return np.array([])

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

    @abstractmethod
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

    @abstractmethod
    def is_stationary(self):
        """Returns whether the kernel is stationary. """


class NormalizedKernelMixin(object):
    """Mixin for kernels which are normalized: k(X, X)=1.

    .. versionadded:: 0.18
    """

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


class StationaryKernelMixin(object):
    """Mixin for kernels which are stationary: k(X, Y)= f(X-Y).

    .. versionadded:: 0.18
    """

    def is_stationary(self):
        """Returns whether the kernel is stationary. """
        return True


class CompoundKernel(Kernel):
    """Kernel which is composed of a set of other kernels.

    .. versionadded:: 0.18
    """

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
        return dict(kernels=self.kernels)

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
            kernel.theta = theta[i * k_dims:(i + 1) * k_dims]

    @property
    def bounds(self):
        """Returns the log-transformed bounds on the theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The log-transformed bounds on the kernel's hyperparameters theta
        """
        return np.vstack([kernel.bounds for kernel in self.kernels])

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
        K_diag : array, shape (n_samples_X, n_kernels)
            Diagonal of kernel k(X, X)
        """
        return np.vstack([kernel.diag(X) for kernel in self.kernels]).T


class KernelOperator(Kernel):
    """Base class for all kernel operators.

    .. versionadded:: 0.18
    """

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
        """Returns the log-transformed bounds on the theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The log-transformed bounds on the kernel's hyperparameters theta
        """
        if self.k1.bounds.size == 0:
            return self.k2.bounds
        if self.k2.bounds.size == 0:
            return self.k1.bounds
        return np.vstack((self.k1.bounds, self.k2.bounds))

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

    .. versionadded:: 0.18

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

    .. versionadded:: 0.18

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

    .. versionadded:: 0.18

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
        """Returns the log-transformed bounds on the theta.

        Returns
        -------
        bounds : array, shape (n_dims, 2)
            The log-transformed bounds on the kernel's hyperparameters theta
        """
        return self.kernel.bounds

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


class ConstantKernel(StationaryKernelMixin, Kernel):
    """Constant kernel.

    Can be used as part of a product-kernel where it scales the magnitude of
    the other factor (kernel) or as part of a sum-kernel, where it modifies
    the mean of the Gaussian process.

    k(x_1, x_2) = constant_value for all x_1, x_2

    .. versionadded:: 0.18

    Parameters
    ----------
    constant_value : float, default: 1.0
        The constant value which defines the covariance:
        k(x_1, x_2) = constant_value

    constant_value_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on constant_value

    """
    def __init__(self, constant_value=1.0, constant_value_bounds=(1e-5, 1e5)):
        self.constant_value = constant_value
        self.constant_value_bounds = constant_value_bounds

    @property
    def hyperparameter_constant_value(self):
        return Hyperparameter(
            "constant_value", "numeric", self.constant_value_bounds)

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

        K = self.constant_value * np.ones((X.shape[0], Y.shape[0]))
        if eval_gradient:
            if not self.hyperparameter_constant_value.fixed:
                return (K, self.constant_value
                        * np.ones((X.shape[0], X.shape[0], 1)))
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
        return self.constant_value * np.ones(X.shape[0])

    def __repr__(self):
        return "{0:.3g}**2".format(np.sqrt(self.constant_value))


class WhiteKernel(StationaryKernelMixin, Kernel):
    """White kernel.

    The main use-case of this kernel is as part of a sum-kernel where it
    explains the noise-component of the signal. Tuning its parameter
    corresponds to estimating the noise-level.

    k(x_1, x_2) = noise_level if x_1 == x_2 else 0

    .. versionadded:: 0.18

    Parameters
    ----------
    noise_level : float, default: 1.0
        Parameter controlling the noise level

    noise_level_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on noise_level

    """
    def __init__(self, noise_level=1.0, noise_level_bounds=(1e-5, 1e5)):
        self.noise_level = noise_level
        self.noise_level_bounds = noise_level_bounds

    @property
    def hyperparameter_noise_level(self):
        return Hyperparameter(
            "noise_level", "numeric", self.noise_level_bounds)

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
            K = self.noise_level * np.eye(X.shape[0])
            if eval_gradient:
                if not self.hyperparameter_noise_level.fixed:
                    return (K, self.noise_level
                            * np.eye(X.shape[0])[:, :, np.newaxis])
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
        return self.noise_level * np.ones(X.shape[0])

    def __repr__(self):
        return "{0}(noise_level={1:.3g})".format(self.__class__.__name__,
                                                 self.noise_level)


class RBF(StationaryKernelMixin, NormalizedKernelMixin, Kernel):
    """Radial-basis function kernel (aka squared-exponential kernel).

    The RBF kernel is a stationary kernel. It is also known as the
    "squared exponential" kernel. It is parameterized by a length-scale
    parameter length_scale>0, which can either be a scalar (isotropic variant
    of the kernel) or a vector with the same number of dimensions as the inputs
    X (anisotropic variant of the kernel). The kernel is given by:

    k(x_i, x_j) = exp(-1 / 2 d(x_i / length_scale, x_j / length_scale)^2)

    This kernel is infinitely differentiable, which implies that GPs with this
    kernel as covariance function have mean square derivatives of all orders,
    and are thus very smooth.

    .. versionadded:: 0.18

    Parameters
    -----------
    length_scale : float or array with shape (n_features,), default: 1.0
        The length scale of the kernel. If a float, an isotropic kernel is
        used. If an array, an anisotropic kernel is used where each dimension
        of l defines the length-scale of the respective feature dimension.

    length_scale_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on length_scale

    """
    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5)):
        self.length_scale = length_scale
        self.length_scale_bounds = length_scale_bounds

    @property
    def anisotropic(self):
        return np.iterable(self.length_scale) and len(self.length_scale) > 1

    @property
    def hyperparameter_length_scale(self):
        if self.anisotropic:
            return Hyperparameter("length_scale", "numeric",
                                  self.length_scale_bounds,
                                  len(self.length_scale))
        return Hyperparameter(
            "length_scale", "numeric", self.length_scale_bounds)

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
        length_scale = _check_length_scale(X, self.length_scale)
        if Y is None:
            dists = pdist(X / length_scale, metric='sqeuclidean')
            K = np.exp(-.5 * dists)
            # convert from upper-triangular matrix to square matrix
            K = squareform(K)
            np.fill_diagonal(K, 1)
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X / length_scale, Y / length_scale,
                          metric='sqeuclidean')
            K = np.exp(-.5 * dists)

        if eval_gradient:
            if self.hyperparameter_length_scale.fixed:
                # Hyperparameter l kept fixed
                return K, np.empty((X.shape[0], X.shape[0], 0))
            elif not self.anisotropic or length_scale.shape[0] == 1:
                K_gradient = \
                    (K * squareform(dists))[:, :, np.newaxis]
                return K, K_gradient
            elif self.anisotropic:
                # We need to recompute the pairwise dimension-wise distances
                K_gradient = (X[:, np.newaxis, :] - X[np.newaxis, :, :]) ** 2 \
                    / (length_scale ** 2)
                K_gradient *= K[..., np.newaxis]
                return K, K_gradient
        else:
            return K

    def __repr__(self):
        if self.anisotropic:
            return "{0}(length_scale=[{1}])".format(
                self.__class__.__name__, ", ".join(map("{0:.3g}".format,
                                                   self.length_scale)))
        else:  # isotropic
            return "{0}(length_scale={1:.3g})".format(
                self.__class__.__name__, np.ravel(self.length_scale)[0])


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

    .. versionadded:: 0.18

    Parameters
    -----------
    length_scale : float or array with shape (n_features,), default: 1.0
        The length scale of the kernel. If a float, an isotropic kernel is
        used. If an array, an anisotropic kernel is used where each dimension
        of l defines the length-scale of the respective feature dimension.

    length_scale_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on length_scale

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
    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5),
                 nu=1.5):
        super(Matern, self).__init__(length_scale, length_scale_bounds)
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
        length_scale = _check_length_scale(X, self.length_scale)
        if Y is None:
            dists = pdist(X / length_scale, metric='euclidean')
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X / length_scale, Y / length_scale,
                          metric='euclidean')

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
            if self.hyperparameter_length_scale.fixed:
                # Hyperparameter l kept fixed
                K_gradient = np.empty((X.shape[0], X.shape[0], 0))
                return K, K_gradient

            # We need to recompute the pairwise dimension-wise distances
            if self.anisotropic:
                D = (X[:, np.newaxis, :] - X[np.newaxis, :, :])**2 \
                    / (length_scale ** 2)
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
                K_gradient = 5.0 / 3.0 * D * (tmp + 1) * np.exp(-tmp)
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
            return "{0}(length_scale=[{1}], nu={2:.3g})".format(
                self.__class__.__name__,
                ", ".join(map("{0:.3g}".format, self.length_scale)),
                self.nu)
        else:
            return "{0}(length_scale={1:.3g}, nu={2:.3g})".format(
                self.__class__.__name__, np.ravel(self.length_scale)[0],
                self.nu)


class RationalQuadratic(StationaryKernelMixin, NormalizedKernelMixin, Kernel):
    """Rational Quadratic kernel.

    The RationalQuadratic kernel can be seen as a scale mixture (an infinite
    sum) of RBF kernels with different characteristic length-scales. It is
    parameterized by a length-scale parameter length_scale>0 and a scale
    mixture parameter alpha>0. Only the isotropic variant where length_scale is
    a scalar is supported at the moment. The kernel given by:

    k(x_i, x_j) = (1 + d(x_i, x_j)^2 / (2*alpha * length_scale^2))^-alpha

    .. versionadded:: 0.18

    Parameters
    ----------
    length_scale : float > 0, default: 1.0
        The length scale of the kernel.

    alpha : float > 0, default: 1.0
        Scale mixture parameter

    length_scale_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on length_scale

    alpha_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on alpha

    """
    def __init__(self, length_scale=1.0, alpha=1.0,
                 length_scale_bounds=(1e-5, 1e5), alpha_bounds=(1e-5, 1e5)):
        self.length_scale = length_scale
        self.alpha = alpha
        self.length_scale_bounds = length_scale_bounds
        self.alpha_bounds = alpha_bounds

    @property
    def hyperparameter_length_scale(self):
        return Hyperparameter(
            "length_scale", "numeric", self.length_scale_bounds)

    @property
    def hyperparameter_alpha(self):
        return Hyperparameter("alpha", "numeric", self.alpha_bounds)

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
            tmp = dists / (2 * self.alpha * self.length_scale ** 2)
            base = (1 + tmp)
            K = base ** -self.alpha
            np.fill_diagonal(K, 1)
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X, Y, metric='sqeuclidean')
            K = (1 + dists / (2 * self.alpha * self.length_scale ** 2)) \
                ** -self.alpha

        if eval_gradient:
            # gradient with respect to length_scale
            if not self.hyperparameter_length_scale.fixed:
                length_scale_gradient = \
                    dists * K / (self.length_scale ** 2 * base)
                length_scale_gradient = length_scale_gradient[:, :, np.newaxis]
            else:  # l is kept fixed
                length_scale_gradient = np.empty((K.shape[0], K.shape[1], 0))

            # gradient with respect to alpha
            if not self.hyperparameter_alpha.fixed:
                alpha_gradient = \
                    K * (-self.alpha * np.log(base)
                         + dists / (2 * self.length_scale ** 2 * base))
                alpha_gradient = alpha_gradient[:, :, np.newaxis]
            else:  # alpha is kept fixed
                alpha_gradient = np.empty((K.shape[0], K.shape[1], 0))

            return K, np.dstack((alpha_gradient, length_scale_gradient))
        else:
            return K

    def __repr__(self):
        return "{0}(alpha={1:.3g}, length_scale={2:.3g})".format(
            self.__class__.__name__, self.alpha, self.length_scale)


class ExpSineSquared(StationaryKernelMixin, NormalizedKernelMixin, Kernel):
    """Exp-Sine-Squared kernel.

    The ExpSineSquared kernel allows modeling periodic functions. It is
    parameterized by a length-scale parameter length_scale>0 and a periodicity
    parameter periodicity>0. Only the isotropic variant where l is a scalar is
    supported at the moment. The kernel given by:

    k(x_i, x_j) = exp(-2 sin(\pi / periodicity * d(x_i, x_j)) / length_scale)^2

    .. versionadded:: 0.18

    Parameters
    ----------
    length_scale : float > 0, default: 1.0
        The length scale of the kernel.

    periodicity : float > 0, default: 1.0
        The periodicity of the kernel.

    length_scale_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on length_scale

    periodicity_bounds : pair of floats >= 0, default: (1e-5, 1e5)
        The lower and upper bound on periodicity

    """
    def __init__(self, length_scale=1.0, periodicity=1.0,
                 length_scale_bounds=(1e-5, 1e5),
                 periodicity_bounds=(1e-5, 1e5)):
        self.length_scale = length_scale
        self.periodicity = periodicity
        self.length_scale_bounds = length_scale_bounds
        self.periodicity_bounds = periodicity_bounds

    @property
    def hyperparameter_length_scale(self):
        return Hyperparameter(
            "length_scale", "numeric", self.length_scale_bounds)

    @property
    def hyperparameter_periodicity(self):
        return Hyperparameter(
            "periodicity", "numeric", self.periodicity_bounds)

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
            arg = np.pi * dists / self.periodicity
            sin_of_arg = np.sin(arg)
            K = np.exp(- 2 * (sin_of_arg / self.length_scale) ** 2)
        else:
            if eval_gradient:
                raise ValueError(
                    "Gradient can only be evaluated when Y is None.")
            dists = cdist(X, Y, metric='euclidean')
            K = np.exp(- 2 * (np.sin(np.pi / self.periodicity * dists)
                              / self.length_scale) ** 2)

        if eval_gradient:
            cos_of_arg = np.cos(arg)
            # gradient with respect to length_scale
            if not self.hyperparameter_length_scale.fixed:
                length_scale_gradient = \
                    4 / self.length_scale**2 * sin_of_arg**2 * K
                length_scale_gradient = length_scale_gradient[:, :, np.newaxis]
            else:  # length_scale is kept fixed
                length_scale_gradient = np.empty((K.shape[0], K.shape[1], 0))
            # gradient with respect to p
            if not self.hyperparameter_periodicity.fixed:
                periodicity_gradient = \
                    4 * arg / self.length_scale**2 * cos_of_arg \
                    * sin_of_arg * K
                periodicity_gradient = periodicity_gradient[:, :, np.newaxis]
            else:  # p is kept fixed
                periodicity_gradient = np.empty((K.shape[0], K.shape[1], 0))

            return K, np.dstack((length_scale_gradient, periodicity_gradient))
        else:
            return K

    def __repr__(self):
        return "{0}(length_scale={1:.3g}, periodicity={2:.3g})".format(
            self.__class__.__name__, self.length_scale, self.periodicity)


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

    .. versionadded:: 0.18

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

    @property
    def hyperparameter_sigma_0(self):
        return Hyperparameter("sigma_0", "numeric", self.sigma_0_bounds)

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
          considered to be a hyperparameter and may be optimized. The other
          kernel parameters are set directly at initialization and are kept
          fixed.

    .. versionadded:: 0.18

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
        self.metric = metric
        self.pairwise_kernels_kwargs = pairwise_kernels_kwargs

    @property
    def hyperparameter_gamma(self):
        return Hyperparameter("gamma", "numeric", self.gamma_bounds)

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
        pairwise_kernels_kwargs = self.pairwise_kernels_kwargs
        if self.pairwise_kernels_kwargs is None:
            pairwise_kernels_kwargs = {}

        X = np.atleast_2d(X)
        K = pairwise_kernels(X, Y, metric=self.metric, gamma=self.gamma,
                             filter_params=True,
                             **pairwise_kernels_kwargs)
        if eval_gradient:
            if self.hyperparameter_gamma.fixed:
                return K, np.empty((X.shape[0], X.shape[0], 0))
            else:
                # approximate gradient numerically
                def f(gamma):  # helper function
                    return pairwise_kernels(
                        X, Y, metric=self.metric, gamma=np.exp(gamma),
                        filter_params=True, **pairwise_kernels_kwargs)
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
        return np.apply_along_axis(self, 1, X).ravel()

    def is_stationary(self):
        """Returns whether the kernel is stationary. """
        return self.metric in ["rbf"]

    def __repr__(self):
        return "{0}(gamma={1}, metric={2})".format(
            self.__class__.__name__, self.gamma, self.metric)
