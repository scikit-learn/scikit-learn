
# This is strongly based on george's kernel module
# Author: Jan Hendrik Metzen <vincent.dubourg@gmail.com>
# Licence: BSD 3 clause

import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform


class Kernel(object):

    def _parse_param_space(self, param_space):
        if not hasattr(param_space, "__iter__"):  # fixed hyperparameter
            self.params = np.array([float(param_space)])
            self.has_bounds = False
            return
        param_space = np.atleast_2d(param_space)
        if param_space.shape[1] == 1:  # fixed hyperparameter
            self.params = param_space[:, 0]
            self.has_bounds = False
        elif param_space.shape[1] == 2:  # lower+upper bound for hyperparameter
            self.bounds = param_space
            self.has_bounds = True
            # Use geometric mean of upper and lower boundary as initial
            # hyperparameter value
            assert not np.any(self.l_bound == None)  # XXX: enforce element-wise comparison to None
            assert not np.any(self.u_bound == None)
            self.params = np.array([np.sqrt(self.l_bound * self.u_bound)])
        elif param_space.shape[1] == 3:  # lower bound, initial value, upper bound
            self.params = param_space[:, 1]
            self.bounds = param_space[:, [0, 2]]
            self.has_bounds = True
        else:
            raise Exception()

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
                                 ", ".join(map("{0}".format,
                                               self.params)))



class KernelOperator(Kernel):

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

    def auto_correlation(self, X, eval_gradient=False):
        if eval_gradient:
            K1, K1_gradient = self.k1.auto_correlation(X, eval_gradient=True)
            K2, K2_gradient = self.k2.auto_correlation(X, eval_gradient=True)
            return K1 + K2, np.dstack((K1_gradient, K2_gradient))
        else:
            return self.k1.auto_correlation(X) + self.k2.auto_correlation(X)

    def cross_correlation(self, X1, X2):
        return self.k1.cross_correlation(X1, X2) \
            + self.k2.cross_correlation(X1, X2)

    def __repr__(self):
        return "{0} + {1}".format(self.k1, self.k2)


class Product(KernelOperator):

    def auto_correlation(self, X, eval_gradient=False):
        if eval_gradient:
            K1, K1_gradient = self.k1.auto_correlation(X, eval_gradient=True)
            K2, K2_gradient = self.k2.auto_correlation(X, eval_gradient=True)
            return K1 * K2, np.dstack((K1_gradient * K2[:, :, None],
                                       K2_gradient * K1[:, :, None]))
        else:
            return self.k1.auto_correlation(X) * self.k2.auto_correlation(X)

    def cross_correlation(self, X1, X2):
        return self.k1.cross_correlation(X1, X2) \
            * self.k2.cross_correlation(X1, X2)

    def __repr__(self):
        return "{0} * {1}".format(self.k1, self.k2)


class ConstantKernel(Kernel):

    def __init__(self, param_space=1.0):
        self._parse_param_space(param_space)

    @property
    def params(self):
        return np.array([self.value])

    @params.setter
    def params(self, theta):
        assert len(theta) == 1
        self.value = theta[0]

    def auto_correlation(self, X, eval_gradient=False):
        K = self.value * np.ones((X.shape[0], X.shape[0]))
        if eval_gradient:
            return K, np.ones((X.shape[0], X.shape[0], 1))
        else:
            return K

    def cross_correlation(self, X1, X2):
        return self.value * np.ones((X1.shape[0], X2.shape[0]))

    def __repr__(self):
        return "{0}".format(self.value)


class RBF(Kernel):
    def __init__(self, param_space=1.0):
        self._parse_param_space(param_space)

    @property
    def params(self):
        return np.asarray(self.l)

    @params.setter
    def params(self, theta):
        self.l = theta

    def auto_correlation(self, X, eval_gradient=False):
        dists = pdist(X / self.l, metric='sqeuclidean')
        K = np.exp(-.5 * dists)
        # convert from upper-triangular matrix to square matrix
        K = squareform(K)
        np.fill_diagonal(K, 1)
        if eval_gradient:
            if self.l.shape[0] == 1:
                K_gradient = (K * squareform(dists) / self.l[0])[:, :, None]
                return K, K_gradient
            elif self.l.shape[0] == X.shape[1]:
                # We need to recompute the pairwise dimension-wise distances
                D = (X[:, None, :] - X[None, :, :])**2 / (self.l ** 3)
                K_gradient = K[..., None] * D
                return K, K_gradient
            else:
                raise Exception("Anisotropic kernels require that the number "
                                "of length scales and features match.")
        else:
            return K

    def cross_correlation(self, X1, X2):
        dists = cdist(X1 / self.l, X2 / self.l, metric='sqeuclidean')
        K = np.exp(-.5 * dists)
        return K


class WhiteKernel(Kernel):
    def __init__(self, param_space=1.0):
        self._parse_param_space(param_space)

    @property
    def params(self):
        return np.asarray([self.c])

    @params.setter
    def params(self, theta):
        self.c = theta[0]

    def auto_correlation(self, X, eval_gradient=False):
        K = self.c * np.eye(X.shape[0])
        if eval_gradient:
            return K, np.eye(X.shape[0])[:, :, None]
        else:
            return K

    def cross_correlation(self, X1, X2):
        K = np.zeros((X1.shape[0], X2.shape[0]))
        K[cdist(X1, X2) < 1e-10] = 1  # entries which are sufficiently similar
        return K
