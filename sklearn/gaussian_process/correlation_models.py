# -*- coding: utf-8 -*-

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# Licence: BSD 3 clause

"""
The built-in correlation models submodule for the gaussian_process module.
"""

import numpy as np

from ..utils import array2d
from ..metrics.pairwise import manhattan_distances


def l1_cross_distances(X):
    """
    Computes the nonzero componentwise L1 cross-distances between the vectors
    in X.

    Parameters
    ----------

    X: array_like
        An array with shape (n_samples, n_features)

    Returns
    -------

    D: array with shape (n_samples * (n_samples - 1) / 2, n_features)
        The array of componentwise L1 cross-distances.

    ij: arrays with shape (n_samples * (n_samples - 1) / 2, 2)
        The indices i and j of the vectors in X associated to the cross-
        distances in D: D[k] = np.abs(X[ij[k, 0]] - Y[ij[k, 1]]).
    """
    X = array2d(X)
    n_samples, n_features = X.shape
    n_nonzero_cross_dist = n_samples * (n_samples - 1) // 2
    ij = np.zeros((n_nonzero_cross_dist, 2), dtype=np.int)
    D = np.zeros((n_nonzero_cross_dist, n_features))
    ll_1 = 0
    for k in range(n_samples - 1):
        ll_0 = ll_1
        ll_1 = ll_0 + n_samples - k - 1
        ij[ll_0:ll_1, 0] = k
        ij[ll_0:ll_1, 1] = np.arange(k + 1, n_samples)
        D[ll_0:ll_1] = np.abs(X[k] - X[(k + 1):n_samples])

    return D, ij


class StationaryCorrelation(object):

    def __init__(self, X, nugget, storage_mode='full'):  # TODO: Storage mode
        self.X = X
        self.nugget = nugget
        self.n_samples = X.shape[0]

        # Calculate array with shape (n_eval, n_features) giving the
        # componentwise distances between locations x and x' at which the
        # correlation model should be evaluated.
        self.D, self.ij = l1_cross_distances(self.X)
        self.D = np.abs(np.asarray(self.D, dtype=np.float))
        # TODO:
#           if (np.min(np.sum(self.D, axis=1)) == 0.
#                    and self.corr != correlation.pure_nugget):
#                raise Exception("Multiple input features cannot have the same"
#                                " value.")


class AbsoluteExponential(StationaryCorrelation):

    def __call__(self, theta, X=None):
        """
        Absolute exponential autocorrelation model.
        (Ornstein-Uhlenbeck stochastic process)::

                                                n
            theta, dx --> r(theta, dx) = exp(  sum  - theta_i * |dx_i| )
                                              i = 1

        Parameters
        ----------
        theta : array_like
            An array with shape 1 (isotropic) or n (anisotropic) giving the
            autocorrelation parameter(s).

        Returns
        -------
        r : array_like
            An array with shape (n_eval, ) containing the values of the
            autocorrelation model.
        """
        theta = np.asarray(theta, dtype=np.float)
        if X is not None:
            # Get pairwise componentwise L1-distances to the input training set
            d = manhattan_distances(X, Y=self.X, sum_over_features=False)
        else:
            d = self.D

        if d.ndim > 1:
            n_features = d.shape[1]
        else:
            n_features = 1

        if theta.size == 1:
            r = np.exp(- theta[0] * np.sum(d, axis=1))
        elif theta.size != n_features:
            raise ValueError("Length of theta must be 1 or %s" % n_features)
        else:
            r = np.exp(- np.sum(theta.reshape(1, n_features) * d, axis=1))

        if X is not None:
            return r.reshape(-1, self.n_samples)
        else:
            R = np.eye(self.n_samples) * (1. + self.nugget)
            R[self.ij[:, 0], self.ij[:, 1]] = r
            R[self.ij[:, 1], self.ij[:, 0]] = r
            return R


class SquaredExponential(StationaryCorrelation):

    def __call__(self, theta, X=None):
        """
        Squared exponential correlation model (Radial Basis Function).
        (Infinitely differentiable stochastic process, very smooth)::

                                                n
            theta, dx --> r(theta, dx) = exp(  sum  - theta_i * (dx_i)^2 )
                                              i = 1

        Parameters
        ----------
        theta : array_like
            An array with shape 1 (isotropic) or n (anisotropic) giving the
            autocorrelation parameter(s).

        Returns
        -------
        r : array_like
            An array with shape (n_eval, ) containing the values of the
            autocorrelation model.
        """
        theta = np.asarray(theta, dtype=np.float)
        if X is not None:
            # Get pairwise componentwise L1-distances to the input training set
            d = manhattan_distances(X, Y=self.X, sum_over_features=False)
        else:
            d = self.D

        if d.ndim > 1:
            n_features = d.shape[1]
        else:
            n_features = 1

        if theta.size == 1:
            r = np.exp(-theta[0] * np.sum(d ** 2, axis=1))
        elif theta.size != n_features:
            raise ValueError("Length of theta must be 1 or %s" % n_features)
        else:
            r = np.exp(-np.sum(theta.reshape(1, n_features) * d ** 2, axis=1))

        if X is not None:
            return r.reshape(-1, self.n_samples)
        else:
            R = np.eye(self.n_samples) * (1. + self.nugget)
            R[self.ij[:, 0], self.ij[:, 1]] = r
            R[self.ij[:, 1], self.ij[:, 0]] = r
            return R


class GeneralizedExponential(StationaryCorrelation):

    def __call__(self, theta, X=None):
        """
        Generalized exponential correlation model.
        (Useful when one does not know the smoothness of the function to be
        predicted.)::

                                                n
            theta, dx --> r(theta, dx) = exp(  sum  - theta_i * |dx_i|^p )
                                              i = 1

        Parameters
        ----------
        theta : array_like
            An array with shape 1+1 (isotropic) or n+1 (anisotropic) giving the
            autocorrelation parameter(s) (theta, p).

        Returns
        -------
        r : array_like
            An array with shape (n_eval, ) with the values of the
            autocorrelation model.
        """
        theta = np.asarray(theta, dtype=np.float)
        if X is not None:
            # Get pairwise componentwise L1-distances to the input training set
            d = manhattan_distances(X, Y=self.X, sum_over_features=False)
        else:
            d = self.D

        if d.ndim > 1:
            n_features = d.shape[1]
        else:
            n_features = 1

        lth = theta.size
        if n_features > 1 and lth == 2:
            theta = np.hstack([np.repeat(theta[0], n_features), theta[1]])
        elif lth != n_features + 1:
            raise Exception("Length of theta must be 2 or %s"
                            % (n_features + 1))
        else:
            theta = theta.reshape(1, lth)

        td = theta[:, 0:-1].reshape(1, n_features) \
            * np.abs(d) ** theta[:, -1]
        r = np.exp(- np.sum(td, 1))

        if X is not None:
            return r.reshape(-1, self.n_samples)
        else:
            R = np.eye(self.n_samples) * (1. + self.nugget)
            R[self.ij[:, 0], self.ij[:, 1]] = r
            R[self.ij[:, 1], self.ij[:, 0]] = r
            return R


class PureNugget(StationaryCorrelation):

    def __call__(self, theta, X=None):
        """
        Spatial independence correlation model (pure nugget).
        (Useful when one wants to solve an ordinary least squares problem!)::

                                                 n
            theta, dx --> r(theta, dx) = 1 if   sum |dx_i| == 0
                                               i = 1
                                         0 otherwise

        Parameters
        ----------
        theta : array_like
            None.


        Returns
        -------
        r : array_like
            An array with shape (n_eval, ) with the values of the
            autocorrelation model.
        """

        theta = np.asarray(theta, dtype=np.float)
        if X is not None:
            # Get pairwise componentwise L1-distances to the input training set
            d = manhattan_distances(X, Y=self.X, sum_over_features=False)
        else:
            d = self.D

        n_eval = d.shape[0]
        r = np.zeros(n_eval)
        r[np.all(d == 0., axis=1)] = 1.

        if X is not None:
            return r.reshape(-1, self.n_samples)
        else:
            R = np.eye(self.n_samples) * (1. + self.nugget)
            R[self.ij[:, 0], self.ij[:, 1]] = r
            R[self.ij[:, 1], self.ij[:, 0]] = r
            return R


class Cubic(StationaryCorrelation):

    def __call__(self, theta, X=None):
        """
        Cubic correlation model::

            theta, dx --> r(theta, dx) =
              n
            prod max(0, 1 - 3(theta_j*d_ij)^2 + 2(theta_j*d_ij)^3) ,  i = 1,...,m
            j = 1

        Parameters
        ----------
        theta : array_like
            An array with shape 1 (isotropic) or n (anisotropic) giving the
            autocorrelation parameter(s).

        Returns
        -------
        r : array_like
            An array with shape (n_eval, ) with the values of the
            autocorrelation model.
        """
        theta = np.asarray(theta, dtype=np.float)
        if X is not None:
            # Get pairwise componentwise L1-distances to the input training set
            d = manhattan_distances(X, Y=self.X, sum_over_features=False)
        else:
            d = self.D

        if d.ndim > 1:
            n_features = self.D.shape[1]
        else:
            n_features = 1

        lth = theta.size
        if lth == 1:
            td = np.abs(d) * theta
        elif lth != n_features:
            raise Exception("Length of theta must be 1 or " + str(n_features))
        else:
            td = np.abs(d) * theta.reshape(1, n_features)

        td[td > 1.] = 1.
        ss = 1. - td ** 2. * (3. - 2. * td)
        r = np.prod(ss, 1)

        if X is not None:
            return r.reshape(-1, self.n_samples)
        else:
            R = np.eye(self.n_samples) * (1. + self.nugget)
            R[self.ij[:, 0], self.ij[:, 1]] = r
            R[self.ij[:, 1], self.ij[:, 0]] = r
            return R


class Linear(StationaryCorrelation):

    def __call__(self, theta, X=None):
        """
        Linear correlation model::

            theta, dx --> r(theta, dx) =
                  n
                prod max(0, 1 - theta_j*d_ij) ,  i = 1,...,m
                j = 1

        Parameters
        ----------
        theta : array_like
            An array with shape 1 (isotropic) or n (anisotropic) giving the
            autocorrelation parameter(s).

        Returns
        -------
        r : array_like
            An array with shape (n_eval, ) with the values of the
            autocorrelation model.
        """
        theta = np.asarray(theta, dtype=np.float)
        if X is not None:
            # Get pairwise componentwise L1-distances to the input training set
            d = manhattan_distances(X, Y=self.X, sum_over_features=False)
        else:
            d = self.D

        if d.ndim > 1:
            n_features = d.shape[1]
        else:
            n_features = 1

        lth = theta.size
        if lth == 1:
            td = np.abs(d) * theta
        elif lth != n_features:
            raise Exception("Length of theta must be 1 or %s" % n_features)
        else:
            td = np.abs(d) * theta.reshape(1, n_features)

        td[td > 1.] = 1.
        ss = 1. - td
        r = np.prod(ss, 1)

        if X is not None:
            return r.reshape(-1, self.n_samples)
        else:
            R = np.eye(self.n_samples) * (1. + self.nugget)
            R[self.ij[:, 0], self.ij[:, 1]] = r
            R[self.ij[:, 1], self.ij[:, 0]] = r
            return R
