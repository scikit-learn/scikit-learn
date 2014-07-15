# -*- coding: utf-8 -*-

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         added factor analysis distance and Matern kernels
# Licence: BSD 3 clause

"""
The built-in correlation models submodule for the gaussian_process module.
"""


import numpy as np


def absolute_exponential(theta, d):
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

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) containing the values of the
        autocorrelation model.
    """
    theta = np.asarray(theta, dtype=np.float)
    d = np.abs(np.asarray(d, dtype=np.float))

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1

    if theta.size == 1:
        return np.exp(- theta[0] * np.sum(d, axis=1))
    elif theta.size != n_features:
        raise ValueError("Length of theta must be 1 or %s" % n_features)
    else:
        return np.exp(- np.sum(theta.reshape(1, n_features) * d, axis=1))


def squared_exponential(theta, d):
    """
    Squared exponential correlation model (Radial Basis Function).
    (Infinitely differentiable stochastic process, very smooth)::

        theta, dx --> r(theta, dx) = exp(-activ),
    where activ=dx.T * M * dx and M is a covariance matrix of size n*n.
    The hyperparameters theta specify
     * a isotropic covariance matrix, i.e., M = theta * I with I being the
       identity, if theta has shape 1
     * an automatic relevance determination model if theta has shape n,
       in which the characteristic length scales of each dimension are learned
       separately:  M = diag(theta)
     * an factor analysis distance model if theta has shape k*n for k> 1,
       in which a low-rank approximation of the full matrix M is learned. This
       low-rank approximation approximates the covariance matrix as low-rank
       matrix plus a diagonal matrix:  M = Lambda * Lambda.T + diag(l),
       where Lambda is a n*(k-1) matrix and l specifies the diagonal matrix.

    See Rasmussen and Williams 2006, p107 for details regarding the different
    variants of the squared exponential kernel.

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic), n (anisotropic), or k*n (factor
        analysis distance) giving the autocorrelation parameter(s). In the
        case of the factor analysis distance, M is approximated by
        M = Lambda * Lambda.T + diag(l), where l is encoded in the last n
        entries of theta and Lambda is encoded row-wise in the first entries of
        theta. Note that Lambda may contain negative entries while theta is
        strictly positive; because of this, the entries of Lambda are set to
        the logarithm with basis 10 of the corresponding entries in theta.

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) containing the values of the
        autocorrelation model.
    """
    return np.exp(-_activation(theta, d))


def matern_1_5(theta, d):
    """
    Matern correlation model vor nu=1.5. Sample paths are once differentiable.

        r(theta, dx) = (1 + np.sqrt(3*activ))*exp(-np.sqrt(3*activ))
    where activ=dx.T * M * dx and M is a covariance matrix of size n*n.
    The hyperparameters theta specify
     * a isotropic covariance matrix, i.e., M = theta * I with I being the
       identity, if theta has shape 1
     * an automatic relevance determination model if theta has shape n,
       in which the characteristic length scales of each dimension are learned
       separately:  M = diag(theta)
     * an factor analysis distance model if theta has shape k*n for k> 1,
       in which a low-rank approximation of the full matrix M is learned. This
       low-rank approximation approximates the covariance matrix as low-rank
       matrix plus a diagonal matrix:  M = Lambda * Lambda.T + diag(l),
       where Lambda is a n*(k-1) matrix and l specifies the diagonal matrix.

    See Rasmussen and Williams 2006, pp84 for details regarding the different
    variants of the Matern kernel.

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic), n (anisotropic), or k*n (factor
        analysis distance) giving the autocorrelation parameter(s). In the
        case of the factor analysis distance, M is approximated by
        M = Lambda * Lambda.T + diag(l), where l is encoded in the last n
        entries of theta and Lambda is encoded row-wise in the first entries of
        theta. Note that Lambda may contain negative entries while theta is
        strictly positive; because of this, the entries of Lambda are set to
        the logarithm with basis 10 of the corresponding entries in theta.

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) containing the values of the
        autocorrelation model.
    """
    activ = _activation(theta, d)
    tmp = np.sqrt(3 * activ)  # temporary variable for preventing recomputation
    return (1 + tmp) * np.exp(-tmp)


def matern_2_5(theta, d):
    """
    Matern correlation model vor nu=2.5. Sample paths are twice differentiable.

       r(theta, dx) = (1 + np.sqrt(5*activ) + 5/3*activ)*exp(-np.sqrt(5*activ))
    where activ=dx.T * M * dx and M is a covariance matrix of size n*n.
    The hyperparameters theta specify
     * a isotropic covariance matrix, i.e., M = theta * I with I being the
       identity, if theta has shape 1
     * an automatic relevance determination model if theta has shape n,
       in which the characteristic length scales of each dimension are learned
       separately:  M = diag(theta)
     * an factor analysis distance model if theta has shape k*n for k> 1,
       in which a low-rank approximation of the full matrix M is learned. This
       low-rank approximation approximates the covariance matrix as low-rank
       matrix plus a diagonal matrix:  M = Lambda * Lambda.T + diag(l),
       where Lambda is a n*(k-1) matrix and l specifies the diagonal matrix.

    See Rasmussen and Williams 2006, pp84 for details regarding the different
    variants of the Matern kernel.

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic), n (anisotropic), or k*n (factor
        analysis distance) giving the autocorrelation parameter(s). In the
        case of the factor analysis distance, M is approximated by
        M = Lambda * Lambda.T + diag(l), where l is encoded in the last n
        entries of theta and Lambda is encoded row-wise in the first entries of
        theta. Note that Lambda may contain negative entries while theta is
        strictly positive; because of this, the entries of Lambda are set to
        the logarithm with basis 10 of the corresponding entries in theta.

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) containing the values of the
        autocorrelation model.
    """
    activ = _activation(theta, d)
    tmp = np.sqrt(5 * activ)  # temporary variable for preventing recomputation
    return (1 + tmp + 5.0 / 3.0 * activ) * np.exp(-tmp)


def generalized_exponential(theta, d):
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

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asarray(theta, dtype=np.float)
    d = np.abs(np.asarray(d, dtype=np.float))

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1

    lth = theta.size
    if n_features > 1 and lth == 2:
        theta = np.hstack([np.repeat(theta[0], n_features), theta[1]])
    elif lth != n_features + 1:
        raise Exception("Length of theta must be 2 or %s" % (n_features + 1))
    else:
        theta = theta.reshape(1, lth)

    td = theta[:, 0:-1].reshape(1, n_features) * d ** theta[:, -1]
    r = np.exp(- np.sum(td, 1))

    return r


def pure_nugget(theta, d):
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

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asarray(theta, dtype=np.float)
    d = np.asarray(d, dtype=np.float)

    n_eval = d.shape[0]
    r = np.zeros(n_eval)
    r[np.all(d == 0., axis=1)] = 1.

    return r


def cubic(theta, d):
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

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asarray(theta, dtype=np.float)
    d = np.abs(np.asarray(d, dtype=np.float))

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1

    lth = theta.size
    if lth == 1:
        td = d * theta
    elif lth != n_features:
        raise Exception("Length of theta must be 1 or " + str(n_features))
    else:
        td = d * theta.reshape(1, n_features)

    td[td > 1.] = 1.
    ss = 1. - td ** 2. * (3. - 2. * td)
    r = np.prod(ss, 1)

    return r


def linear(theta, d):
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

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the values of the autocorrelation
        model.
    """

    theta = np.asarray(theta, dtype=np.float)
    d = np.abs(np.asarray(d, dtype=np.float))

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1

    lth = theta.size
    if lth == 1:
        td = d * theta
    elif lth != n_features:
        raise Exception("Length of theta must be 1 or %s" % n_features)
    else:
        td = d * theta.reshape(1, n_features)

    td[td > 1.] = 1.
    ss = 1. - td
    r = np.prod(ss, 1)

    return r


def _activation(theta, dx):
    """ Utility function for computing activation in correlation models.

    Computes the activation activ=dx.T * M * dx where M is a covariance matrix
    of size n*n. The hyperparameters theta specify
     * an isotropic covariance matrix, i.e., M = theta * I with I being the
       identity, if theta has shape 1
     * an automatic relevance determination model if theta has shape n,
       in which the characteristic length scales of each dimension are learned
       separately:  M = diag(theta)
     * a factor analysis distance model if theta has shape k*n for k> 1,
       in which a low-rank approximation of the full matrix M is learned. This
       low-rank approximation approximates the covariance matrix as low-rank
       matrix plus a diagonal matrix:  M = Lambda * Lambda.T + diag(l),
       where Lambda is a n*(k-1) matrix and l specifies the diagonal matrix.

    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic), n (anisotropic), or k*n (factor
        analysis distance) giving the autocorrelation parameter(s). In the
        case of the factor analysis distance, M is approximated by
        M = Lambda * Lambda.T + diag(l), where l is encoded in the last n
        entries of theta and Lambda is encoded row-wise in the first entries of
        theta. Note that Lambda may contain negative entries while theta is
        strictly positive; because of this, the entries of Lambda are set to
        the logarithm with basis 10 of the corresponding entries in theta.

    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        differences of x and x' at which the correlation model
        should be evaluated.

    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) with the activation values for the
        respective componentwise differences dx.
    """
    theta = np.asarray(theta, dtype=np.float)
    dx = np.asarray(dx, dtype=np.float)

    if dx.ndim > 1:
        n_features = dx.shape[1]
    else:
        n_features = 1

    if theta.size == 1:  # case where M is isotropic: M = diag(theta[0])
        return theta[0] * np.sum(dx ** 2, axis=1)
    elif theta.size == n_features:  # anisotropic but diagonal case (ARD)
        return np.sum(theta.reshape(1, n_features) * dx ** 2, axis=1)
    elif theta.size % n_features == 0:
        # Factor analysis case: M = lambda*lambda.T + diag(l)
        theta = theta.reshape((1, theta.size))
        M = np.diag(theta[0, :n_features])  # the diagonal matrix part l
        # The low-rank matrix contribution which allows accounting for
        # correlations in the feature dimensions
        # NOTE: these components of theta are passed through a log-function to
        #       allow negative values in Lambda
        Lambda = np.log10(theta[0, n_features:].reshape((n_features, -1)))
        M += Lambda.dot(Lambda.T)
        return np.sum(dx.dot(M) * dx, -1)
    else:
        raise ValueError("Length of theta must be 1 or a multiple of %s."
                         % n_features)
