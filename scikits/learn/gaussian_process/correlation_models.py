#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# License: BSD style

"""
The built-in regression models submodule for the gaussian_process module.
"""

################
# Dependencies #
################

import numpy as np


#############################
# Defaut correlation models #
#############################


def absolute_exponential(theta, d):
    """
    Absolute exponential autocorrelation model.
    (Ornstein-Uhlenbeck stochastic process)

                                        n
    theta, dx --> r(theta, dx) = exp(  sum  - theta_i * |dx_i| )
                                      i = 1

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic) or n (anisotropic) giving the
        autocorrelation parameter(s).

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        distances between locations x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) containing the values of the
        autocorrelation model.
    """

    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1
    if theta.size == 1:
        theta = np.repeat(theta, n_features)
    elif theta.size != n_features:
        raise ValueError("Length of theta must be 1 or " + str(n_features))

    td = - theta.reshape(1, n_features) * abs(d)
    r = np.exp(np.sum(td, 1))

    return r


def squared_exponential(theta, d):
    """
    Squared exponential correlation model (Radial Basis Function).
    (Infinitely differentiable stochastic process, very smooth)

                                        n
    theta, dx --> r(theta, dx) = exp(  sum  - theta_i * (dx_i)^2 )
                                      i = 1

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic) or n (anisotropic) giving the
        autocorrelation parameter(s).

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        distances between locations x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) containing the values of the
        autocorrelation model.
    """

    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1
    if theta.size == 1:
        theta = np.repeat(theta, n_features)
    elif theta.size != n_features:
        raise Exception("Length of theta must be 1 or " + str(n_features))

    td = - theta.reshape(1, n_features) * d ** 2
    r = np.exp(np.sum(td, 1))

    return r


def generalized_exponential(theta, d):
    """
    Generalized exponential correlation model.
    (Useful when one does not know the smoothness of the function to be
    predicted.)

                                        n
    theta, dx --> r(theta, dx) = exp(  sum  - theta_i * |dx_i|^p )
                                      i = 1

    Parameters
    ----------
    theta : array_like
        An array with shape 1+1 (isotropic) or n+1 (anisotropic) giving the
        autocorrelation parameter(s) (theta, p).

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        distances between locations x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1
    lth = theta.size
    if n_features > 1 and lth == 2:
        theta = np.hstack([np.repeat(theta[0], n_features), theta[1]])
    elif lth != n_features + 1:
        raise Exception("Length of theta must be 2 or " + str(n_features + 1))
    else:
        theta = theta.reshape(1, lth)

    td = - theta[:, 0:-1].reshape(1, n_features) * abs(d) ** theta[:, -1]
    r = np.exp(np.sum(td, 1))

    return r


def pure_nugget(theta, d):
    """
    Spatial independence correlation model (pure nugget).
    (Useful when one wants to solve an ordinary least squares problem!)

                                         n
    theta, dx --> r(theta, dx) = 1 if   sum |dx_i| == 0
                                       i = 1
                                 0 otherwise

    Parameters
    ----------
    theta : array_like
        None.

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        distances between locations x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)

    n_eval = d.shape[0]
    r = np.zeros(n_eval)
    # The ones on the diagonal of the correlation matrix are enforced within
    # the KrigingModel instanciation to allow multiple design sites in this
    # ordinary least squares context.

    return r


def cubic(theta, d):
    """
    Cubic correlation model.

    theta, dx --> r(theta, dx) =
          n
        prod max(0, 1 - 3(theta_j*d_ij)^2 + 2(theta_j*d_ij)^3) ,  i = 1,...,m
        j = 1

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic) or n (anisotropic) giving the
        autocorrelation parameter(s).

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        distances between locations x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1
    lth = theta.size
    if  lth == 1:
        theta = np.repeat(theta, n_features)[np.newaxis][:]
    elif lth != n_features:
        raise Exception("Length of theta must be 1 or " + str(n_features))
    else:
        theta = theta.reshape(1, n_features)

    td = abs(d) * theta
    td[td > 1.] = 1.
    ss = 1. - td ** 2. * (3. - 2. * td)
    r = np.prod(ss, 1)

    return r


def linear(theta, d):
    """
    Linear correlation model.

    theta, dx --> r(theta, dx) =
          n
        prod max(0, 1 - theta_j*d_ij) ,  i = 1,...,m
        j = 1

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic) or n (anisotropic) giving the
        autocorrelation parameter(s).

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        distances between locations x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1
    lth = theta.size
    if  lth == 1:
        theta = np.repeat(theta, n_features)[np.newaxis][:]
    elif lth != n_features:
        raise Exception("Length of theta must be 1 or " + str(n_features))
    else:
        theta = theta.reshape(1, n_features)

    td = abs(d) * theta
    td[td > 1.] = 1.
    ss = 1. - td
    r = np.prod(ss, 1)

    return r
