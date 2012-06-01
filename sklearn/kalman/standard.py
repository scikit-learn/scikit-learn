"""
=====================================
Inference for Linear-Gaussian Systems
=====================================

This module implements the Kalman Filter, Kalman Smoother, and
EM Algorithm for Linear-Gaussian state space models.  In other
words, assuming that

.. math::

    x_{t+1}   &= A_t x_t + b_t + e_t^1        \\
    z_{t}     &= C_t x_t + d_t + e_t^2        \\
    e_t^1     &\sim \text{MultivariateNormal}(0, Q)       \\
    e_t^2     &\sim \text{MultivariateNormal}(0, R)       \\
    x_0       &\sim \text{MultivariateNormal}(\mu_0, \Sigma_0)

then the Kalman Filter calculates exactly,

.. math::

    P(x_t | z_{0:t})

the Kalman Smoother,

.. math::

    P(x_t | z_{0:T-1})

the EM algorithm, for :math:`\theta = (A,b,C,d,Q,R,\mu_0,\Sigma_0)`,

.. math::

    \arg\max_{\theta} P(z_{0:T-1}; \theta)

Observations and states are assumed to be sampled at times [0...T-1]

References
==========
    * Abbeel, Pieter. "Maximum Likelihood, EM".
      http://www.cs.berkeley.edu/~pabbeel/cs287-fa11/
    * Yu, Byron M. and Shenoy, Krishna V. and Sahani, Maneesh. "Derivation of
      Kalman Filtering and Smoothing Equations".
      http://www.ece.cmu.edu/~byronyu/papers/derive_ks.pdf
    * Ghahramani, Zoubin and Hinton, Geoffrey E. "Parameter Estimation for Linear
      Dynamical Systems." http://mlg.eng.cam.ac.uk/zoubin/course04/tr-96-2.pdf
    * Welling, Max. "The Kalman Filter".
      http://www.cs.toronto.edu/~welling/classnotes/papers_class/KF.ps.gz
"""
import warnings

import numpy as np
import numpy.ma as ma

from ..base import BaseEstimator
from ..mixture import log_multivariate_normal_density
from ..utils import array1d, array2d

# Dimensionality of each Kalman Filter parameter for a single time step
DIM = {
    'A': 2,
    'b': 1,
    'C': 2,
    'd': 1,
    'Q': 2,
    'R': 2,
    'mu_0': 1,
    'sigma_0': 2,
}

def _last_dims(X, t, ndims=2):
    """Extract the final `ndim` dimensions at index `t` if
    X has >= ndim+1 dimensions, else return X.

    Parameters
    ==========
    X : array with at least dimension `ndims`
    t : int
        index to use for the `ndims`+1th dimension
    ndims : int, optional
        number of dimensions in the array desired

    Returns
    =======
    Y : array with dimension `ndims`
        the final `ndims` dimensions indexed by `t`
    """
    X = np.asarray(X)
    if len(X.shape) == ndims + 1:
        return X[t]
    elif len(X.shape) == ndims:
        return X
    else:
        raise ValueError(("X only has %d dimensions when %d" + \
                " or more are required") % (len(X.shape), ndims))


def _pad(X, n=1, dim=None):
    """
    Pad an array along its first axis with `n` all-zero sub-arrays if its
    dimension isn't equal to `dim`.
    """
    if len(X.shape) == dim:
        return X
    else:
        xs = [np.zeros( X.shape[1:] )[np.newaxis] for i in range(n)]
        return np.vstack(xs + [X])


def _logmvnpdf(x, mu, sigma):
    """log density of the multivariate normal distribution

    Parameters
    ==========
    x : [n_dim] array
        point sampled from multivariate normal
    mu : [n_dim] array
        mean of multivariate normal
    sigma : [n_dim, n_dim] array
        covariance of multivariate normal
    """
    return log_multivariate_normal_density(x[np.newaxis,:], mu[np.newaxis,:],
        sigma[np.newaxis,:,:], 'full')


def _filter_predict(A, Q, b, mu, sigma):
    """Calculate the mean and covariance of :math:`P(x_{t+1} | z_{0:t})`

    Parameters
    ==========
    A : [n_dim_state, n_dim_state} array
        state transition matrix from time t to t+1
    Q : [n_dim_state, n_dim_state] array
        covariance matrix for state transition from time t to t+1
    b : [n_dim_state] array
        offset for state transition from time t to t+1
    mu: [n_dim_state] array
        mean of state at time t given observations from times
        [0...t]
    sigma: [n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [0...t]

    Returns
    =======
    mu_pred : [n_dim_state] array
        mean of state at time t+1 given observations from times [0...t]
    sigma_pred : [n_dim_state, n_dim_state] array
        covariance of state at time t+1 given observations from times
        [0...t]
    """
    mu_pred = A.dot(mu) + b
    sigma_pred = A.dot(sigma).dot(A.T) + Q

    return (mu_pred, sigma_pred)


def _filter_correct(C, R, d, mu_pred, sigma_pred, z):
    """Incorporate observation z from time t to turn :math:`P(x_t | z_{0:t-1})`
    into :math:`P(x_t | z_{0:t})`

    Parameters
    ==========
    C : [n_dim_obs, n_dim_state] array
        observation matrix for time t
    R : [n_dim_obs, n_dim_obs] array
        covariance matrix for observation at time t
    d : [n_dim_obs] array
        offset for observation at time t
    mu_pred : [n_dim_state] array
        mean of state at time t given observations from times
        [0...t-1]
    sigma_pred : [n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [0...t-1]
    z : [n_dim_obs] array
        observation at time t.  If `z` is a masked array and any of its
        values are masked, the observation will be ignored.

    Returns
    =======
    K : [n_dim_state, n_dim_obs] array
        Kalman gain matrix for time t
    mu_new : [n_dim_state] array
        mean of state at time t given observations from times
        [0...t]
    sigma_new : [n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [0...t]
    ll : float
        log likelihood of observation at time t given observations
        from times [0...t-1]
    """
    if not np.any(ma.getmask(z)):
        mu_obs = C.dot(mu_pred) + d
        sigma_obs = C.dot(sigma_pred).dot(C.T) + R
        ll = _logmvnpdf(z, mu_obs, sigma_obs)

        K_new = sigma_pred.dot(C.T).dot(np.linalg.pinv(sigma_obs))

        mu_new = mu_pred + K_new.dot(z - mu_obs)
        sigma_new = sigma_pred - K_new.dot(C).dot(sigma_pred)
    else:
        n_dim_state = sigma_pred.shape[0]
        n_dim_obs = C.shape[0]
        K_new = np.zeros( (n_dim_state, n_dim_obs) )

        ll = 0.0
        mu_new = mu_pred
        sigma_new = sigma_pred

    return (K_new, mu_new, sigma_new, ll)


def _filter(A, C, Q, R, b, d, mu_0, sigma_0, Z):
    """Apply the Kalman Filter to calculate posterior distribution over hidden
    states given observations up to and including the current time step.

    Parameters
    ==========
    A : [T-1, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
        state transition matrix
    C : [T, n_dim_obs, n_dim_obs] or [n_dim_obs, n_dim_obs] array-like
        observation matrix
    Q : [T-1, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
        state transition covariance matrix
    R : [T, n_dim_obs, n_dim_obs] or [n_dim_obs, n_dim_obs] array-like
        observation covariance matrix
    b : [T-1, n_dim_state] or [n_dim_state] array-like
        state offset
    d : [T, n_dim_obs] or [n_dim_obs] array-like
        observations for times [0...T-1]
    mu_0 : [n_dim_state] array-like
        mean of initial state distribution
    sigma_0 : [n_dim_state, n_dim_state] array-like
        covariance of initial state distribution
    Z : [T, n_dim_obs] array
        observations from times [0...T-1].  If `Z` is a masked array and any of
        `Z[t]` is masked, then `Z[t]` will be treated as a missing observation.

    Returns
    =======
    mu_pred : [T, n_dim_state] array
        `mu_pred[t]` = mean of hidden state at time t given observations from
        times [0...t-1]
    sigma_pred : [T, n_dim_state, n_dim_state] array
        `sigma_pred[t]` = covariance of hidden state at time t given
        observations from times [0...t-1]
    K : [T, n_dim_state] array
        `K[t]` = Kalman gain matrix for time t
    mu_filt : [T, n_dim_state] array
        `mu_filt[t]` = mean of hidden state at time t given observations from
        times [0...t]
    sigma_filt : [T, n_dim_state] array
        `sigma_filt[t]` = covariance of hidden state at time t given
        observations from times [0...t]
    ll : [T] array
        `ll[t]` = log likelihood of observation at time t given observations
        from times [0...t-1]
    """
    T = Z.shape[0]
    n_dim_state = len(mu_0)
    n_dim_obs = Z.shape[1]

    mu_pred = np.zeros((T, n_dim_state))
    sigma_pred = np.zeros((T, n_dim_state, n_dim_state))
    K = np.zeros((T, n_dim_state, n_dim_obs))
    mu_filt = np.zeros((T, n_dim_state))
    sigma_filt = np.zeros((T, n_dim_state, n_dim_state))
    ll = np.zeros(T)

    for t in range(T):
        if t == 0:
            mu_pred[t] = mu_0
            sigma_pred[t] = sigma_0
        else:
            A_t1 = _last_dims(A, t - 1)
            Q_t1 = _last_dims(Q, t - 1)
            b_t1 = _last_dims(b, t - 1, ndims=1)
            (mu_pred[t], sigma_pred[t]) = _filter_predict(
                A_t1, Q_t1, b_t1, mu_filt[t-1], sigma_filt[t-1])

        C_t = _last_dims(C, t)
        R_t = _last_dims(R, t)
        d_t = _last_dims(d, t, ndims=1)
        (K[t], mu_filt[t], sigma_filt[t], ll[t]) = _filter_correct(
            C_t, R_t, d_t, mu_pred[t], sigma_pred[t], Z[t])

    return (mu_pred, sigma_pred, K, mu_filt, sigma_filt, ll)


def _smooth_update(A, mu_for, sigma_for, mu_pred, sigma_pred, mu_rev, sigma_rev):
    """One iteration of Kalman Smoothing.  Calculates posterior
    distribution of the hidden state at time t given the observations
    from times [1...T].

    Parameters
    ==========
    A : [n_dim_state, n_dim_state] array
        state transition matrix from time t to t+1
    mu_for : [n_dim_state] array
        mean of filtered state at time t given observations from
        times [1...t]
    sigma_for : [n_dim_state, n_dim_state] array
        covariance of filtered state at time t given observations from
        times [1...t]
    mu_pred : [n_dim_state] array
        mean of filtered state at time t+1 given observations from
        times [1...t]
    sigma_pred : [n_dim_state, n_dim_state] array
        covariance of filtered state at time t+1 given observations from
        times [1...t]
    mu_rev : [n_dim_state] array
        mean of smoothed state at time t+1 given observations from
        times [1...T]
    sigma_rev : [n_dim_state, n_dim_state] array
        covariance of smoothed state at time t+1 given observations from
        times [1...T]

    Returns
    =======
    mu : [n_dim_state] array
        mean of smoothed state at time t given observations from times
        [1...T]
    sigma : [n_dim_state, n_dim_state] array
        covariance of smoothed state at time t given observations from
        times [1...T]
    L : [n_dim_state, n_dim_state] array
        correction matrix for Kalman Smoothing at time t
    """
    L = sigma_for.dot(A.T).dot(np.linalg.pinv(sigma_pred))

    mu = mu_for + L.dot(mu_rev - mu_pred)
    sigma = sigma_for + L.dot(sigma_rev - sigma_pred).dot(L.T)

    return (mu, sigma, L)


def _smooth(A, mu_filt, sigma_filt, mu_pred, sigma_pred):
    """Apply the Kalman Smoother to estimate the hidden state at time
    t for t = [0...T-1] given all observations.

    Parameters
    ==========
    A : [T-1, n_dim_state, n_dim_state]  or [n_dim_state, n_dim_state] array
        `A[t]` = transition matrix from time t to t+1
    mu_filt : [T, n_dim_state] array
        `mu_filt[t]` = mean state estimate for time t given observations from
        times [0...t]
    sigma_filt : [T, n_dim_state, n_dim_state] array
        `sigma_filt[t]` = covariance of state estimate for time t given
        observations from times [0...t]
    mu_pred : [T, n_dim_state] array
        `mu_pred[t]` = mean state estimate for time t given observations from
        times [0...t-1]
    sigma_pred : [T, n_dim_state, n_dim_state] array
        `sigma_pred[t]` = covariance of state estimate for time t given
        observations from times [0...t-1]

    Returns
    =======
    mu_smooth : [T, n_dim_state]
        mean of hidden state distributions for times [0...T-1] given
        all observations
    sigma_smooth : [T, n_dim_state, n_dim_state] array
        covariance matrix of hidden state distributions for times
        [0...T-1] given all observations
    L : [T-1, n_dim_state, n_dim_state] array
        Kalman Smoothing correction matrices for times [0...T-2]
    """
    T, n_dim_state = mu_filt.shape

    mu_smooth = np.zeros((T, n_dim_state))
    sigma_smooth = np.zeros((T, n_dim_state, n_dim_state))
    L = np.zeros((T - 1, n_dim_state, n_dim_state))

    mu_smooth[-1] = mu_filt[-1]
    sigma_smooth[-1] = sigma_filt[-1]

    for t in reversed(range(T-1)):
        A_t = _last_dims(A, t)
        (mu_smooth[t], sigma_smooth[t], L[t]) = _smooth_update(
            A_t,
            mu_filt[t], sigma_filt[t],
            mu_pred[t + 1], sigma_pred[t + 1],
            mu_smooth[t + 1], sigma_smooth[t + 1],
        )
    return (mu_smooth, sigma_smooth, L)


def _smooth_pair(sigma_smooth, L):
    """Calculate covariance between hidden states at t and t-1 for t = [0...T-1]

    Parameters
    ==========
    sigma_smooth : [T, n_dim_state, n_dim_state] array
        covariance of hidden state given all observations
    L : [T-1, n_dim_state, n_dim_state]
        Correction matrices from Kalman Smoothing

    Returns
    =======
    sigma_pair_smooth : [T, n_dim_state, n_dim_state] array
        Covariance between hidden states at times t and t-1 for t = [1...T-1].
        Time 0 is ignored.
    """
    T, n_dim_state, _ = sigma_smooth.shape
    sigma_pair_smooth = np.zeros((T, n_dim_state, n_dim_state))
    for t in range(1, T):
        sigma_pair_smooth[t] = sigma_smooth[t].dot(L[t - 1].T)
    return sigma_pair_smooth


def _em(Z, b, d, mu_smooth, sigma_smooth, sigma_pair_smooth, given={}):
    """Estimate Linear-Gaussian model parameters by maximizing the expected log
    likelihood.

    Parameters
    ==========
    Z : [T, n_dim_obs] array
        observations for times [0...T-1].  If Z is a masked array and any of
        Z[t] is masked, then it will be treated as a missing observation.
    b : [n_dim_state] or [T-1, n_dim_state] array
        transition offset
    d : [n_dim_obs] or [T, n_dim_obs] array
        observation offsets
    mu_smooth : [T, n_dim_state] array
        mu_smooth[t] = mean of state at time t given all observations
    sigma_smooth : [T, n_dim_state, n_dim_state] array
        sigma_smooth[t] = covariance of state at time t given all observations
    sigma_pair_smooth : [T, n_dim_state, n_dim_state] array
        sigma_pair_smooth[t] = covariance between states at times t and
        t-1 given all observations.  sigma_pair_smooth[0] is ignored.
    given: dict
        if one of the variables EM is capable of predicting is in given, then
        that value will be used and EM will not attempt to estimate it.  e.g.,
        if 'C' is in given, C will not be estimated and given['C'] will be
        returned in its place.

    Returns
    =======
    A : [n_dim_state, n_dim_state] array
        estimated transition matrix
    C : [n_dim_obs, n_dim_state] array
        estimated observation matrix
    b : [n_dim_state] array
        estimated transition offset
    d : [n_dim_obs] array
        estimated observation offset
    Q : [n_dim_state, n_dim_state] array
        estimated covariance matrix for state transitions
    R : [n_dim_obs, n_dim_obs] array
        estimated covariance matrix for observations
    mu_0 : [n_dim_state] array
        estimated mean of initial state distribution
    sigma_0 : [n_dim_state] array
        estimated covariance of initial state distribution
    """
    C = given.get('C', _em_C(Z, d, mu_smooth, sigma_smooth))
    R = given.get('R', _em_R(Z, d, C, mu_smooth, sigma_smooth))
    A = given.get('A', _em_A(b, mu_smooth, sigma_smooth, sigma_pair_smooth))
    Q = given.get('Q', _em_Q(A, b, mu_smooth, sigma_smooth, sigma_pair_smooth))
    mu_0 = given.get('mu_0', _em_mu_0(mu_smooth))
    sigma_0 = given.get('sigma_0', _em_sigma_0(mu_0, mu_smooth, sigma_smooth))
    b = given.get('b', _em_b(A, mu_smooth))
    d = given.get('d', _em_d(C, mu_smooth, Z))

    return (A, C, b, d, Q, R, mu_0, sigma_0)


def _em_C(Z, d, mu_smooth, sigma_smooth):
    """Maximize expected log likelihood of observations with respect to the
    observation matrix C.

    .. math::

        C &= ( \sum_{t=0}^{T-1} (y_t - d_t) \mathbb{E}[x_t] )
             ( \sum_{t=0}^{T-1} \mathbb{E}[x_t x_t^T] )^-1
    """
    _, n_dim_state = mu_smooth.shape
    T, n_dim_obs = Z.shape
    res1 = np.zeros((n_dim_obs, n_dim_state))
    res2 = np.zeros((n_dim_state, n_dim_state))
    for t in range(T):
        if not np.any(ma.getmask(Z[t])):
            d_t = _last_dims(d, t, ndims=1)
            res1 += np.outer(Z[t] - d_t, mu_smooth[t])
            res2 += sigma_smooth[t] + np.outer(mu_smooth[t], mu_smooth[t])
    return res1.dot(np.linalg.pinv(res2))


def _em_R(Z, d, C, mu_smooth, sigma_smooth):
    """Maximize expected log likelihood of observations with respect to the
    observation covariance matrix R.

    .. math::

        R &= \frac{1}{T} \sum_{t=0}^{T-1}
                [y_t - C_t \mathbb{E}[x_t] - d_t]
                    [y_t - C_t \mathbb{E}[x_t] - d_t]^T
                + C_t Var(x_t) C_t^T
    """
    _, n_dim_state = mu_smooth.shape
    T, n_dim_obs = Z.shape
    res = np.zeros((n_dim_obs, n_dim_obs))
    n_obs = 0
    for t in range(T):
        if not np.any(ma.getmask(Z[t])):
            C_t = _last_dims(C, t)
            d_t = _last_dims(d, t, ndims=1)
            err = Z[t] - C_t.dot(mu_smooth[t]) - d_t
            res += np.outer(err, err) \
                + C_t.dot(sigma_smooth[t]).dot(C_t.T)
            n_obs += 1
    if n_obs > 0:
        return (1.0 / n_obs) * res
    else:
        return res


def _em_A(b, mu_smooth, sigma_smooth, sigma_pair_smooth):
    """Maximize expected log likelihood of observations with respect to the
    transition matrix A.

    .. math::

        A &= ( \sum_{t=1}^{T-1} \mathbb{E}[x_t x_{t-1}^{T}]
                - b_{t-1} \mathbb{E}[x_{t-1}]^T )
             ( \sum_{t=1}^{T-1} \mathbb{E}[x_{t-1} x_{t-1}^T] )^{-1}
    """
    T, n_dim_state, _ = sigma_smooth.shape
    res1 = np.zeros((n_dim_state, n_dim_state))
    res2 = np.zeros((n_dim_state, n_dim_state))
    for t in range(1, T):
        b_t1 = _last_dims(b, t - 1, ndims=1)
        res1 += sigma_pair_smooth[t] \
            + np.outer(mu_smooth[t], mu_smooth[t - 1]) \
            - np.outer(b_t1, mu_smooth[t - 1])
        res2 += sigma_smooth[t - 1] \
            + np.outer(mu_smooth[t - 1], mu_smooth[t - 1])
    return res1.dot(np.linalg.pinv(res2))


def _em_Q(A, b, mu_smooth, sigma_smooth, sigma_pair_smooth):
    """Maximize expected log likelihood of observations with respect to the
    transition covariance matrix Q.

    .. math::

        Q &= \frac{1}{T-1} \sum_{t=0}^{T-2}
                (\mathbb{E}[x_{t+1}] - A_t \mathbb{E}[x_t] - b_t)
                    (\mathbb{E}[x_{t+1}] - A_t \mathbb{E}[x_t] - b_t)^T
                + A_t Var(x_t) A_t^T + Var(x_{t+1})
                - Cov(x_{t+1}, x_t) A_t^T - A_t Cov(x_t, x_{t+1})
    """
    T, n_dim_state, _ = sigma_smooth.shape
    res = np.zeros((n_dim_state, n_dim_state))
    for t in range(T - 1):
        A_t = _last_dims(A, t)
        b_t = _last_dims(b, t, ndims=1)
        err = mu_smooth[t + 1] - A_t.dot(mu_smooth[t]) - b_t
        Vt1t_A = sigma_pair_smooth[t + 1].dot(A_t.T)
        res += np.outer(err, err) + A_t.dot(sigma_smooth[t]).dot(A_t.T) \
            + sigma_smooth[t + 1] - Vt1t_A - Vt1t_A.T

    return (1.0 / (T-1)) * res


def _em_mu_0(mu_smooth):
    """Maximize expected log likelihood of observations with respect to the
    initial state distribution mean mu_0.

    .. math::

        mu_0 = \mathbb{E}[x_0]
    """

    return mu_smooth[0]


def _em_sigma_0(mu_0, mu_smooth, sigma_smooth):
    """Maximize expected log likelihood of observations with respect to the
    covariance of the initial state distribution sigma_0

    .. math::

        \sigma_0 = \mathbb{E}[x_0, x_0^T] - mu_0 \mathbb{E}[x_0]^T
                   - \mathbb{E}[x_0] mu_0^T + mu_0 mu_0^T
    """
    x0 = mu_smooth[0]
    x0_x0 = sigma_smooth[0] + np.outer(x0, x0)
    return x0_x0 - np.outer(mu_0, x0) - np.outer(x0, mu_0) \
        + np.outer(mu_0, mu_0)


def _em_b(A, mu_smooth):
    """Maximize expected log likelihood of observations with respect to the
    state transition offset

    .. math::

        b = \frac{1}{T-1} \sum_{t=1}^{T-1}
                \mathbb{E}[x_t] - A_{t-1} \mathbb{E}[x_{t-1}]
    """
    T, n_dim_state = mu_smooth.shape
    b = np.zeros(n_dim_state)
    for t in range(1, T):
        A_t1 = _last_dims(A, t-1)
        b += mu_smooth[t] - A_t1.dot(mu_smooth[t - 1])
    return (1.0/(T-1))*b


def _em_d(C, mu_smooth, Z):
    """Maximize expected log likelihood of observations with respect to the
    observation offset

    .. math::

        d = \frac{1}{T} \sum_{t=0}^{T-1} z_t - C_{t} \mathbb{E}[x_{t}]
    """
    T, n_dim_obs = Z.shape
    d = np.zeros(n_dim_obs)
    n_obs = 0
    for t in range(T):
        if not np.any(ma.getmask(Z[t])):
            C_t = _last_dims(C, t)
            d += Z[t] - C_t.dot(mu_smooth[t])
            n_obs += 1
    if n_obs > 0:
        return (1.0/n_obs)*d
    else:
        return d


class KalmanFilter(BaseEstimator):
    """
    Implements the Kalman Filter, Kalman Smoother, and EM algorithm for
    a Linear Gaussian model specified by,

    .. math::

        x_{t+1}   &= A_{t} x_{t} + b_{t} + Normal(0, Q_{t}) \\
        z_{t}     &= C_{t} x_{t} + d_{t} + Normal(0, R_{t})

    The Kalman Filter is an algorithm designed to estimate 
    :math:`P(x_t | z_{0:t})``.  As all state transitions and observations are
    linear with Gaussian distributed noise, these distributions can be
    represented exactly as Gaussian distributions with mean `mu_filt[t]` and
    covariances `sigma_filt[t]`.

    Similarly, the Kalman Smoother is an algorithm designed to estimate
    :math:`P(x_t | z_{0:T})`.

    The EM algorithm aims to find for :math:`theta = (A, b, C, d, Q, R, \mu_0,
    \sigma_0)`

    .. math::

        \max_{\theta} P(z_{0:T-1}; \theta)

    If we define :math:`L(x_{0:t},\theta) = \log P(z_{0:T-1}, x_{0:T-1};
    \theta)`, then the EM algorithm works by iteratively finding,

    .. math::

        P(x_{0:T-1} | z_{0:T-1}, \theta_i)

    then by maximizing,

    .. math::

        \theta_{i+1} = \arg\max_{\theta}
            E_{x_{0:T-1}} [ L(x_{0:t}, \theta)| z_{0:T-1}, \theta_i ]

    Parameters
    ==========
    A : [T-1, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
        state transition matrix for times [0...T-2]
    C : [T, n_dim_obs, n_dim_obs] or [n_dim_obs, n_dim_obs] array-like
        observation matrix for times [0...T-1]
    Q : [n_dim_state, n_dim_state] array-like
        state transition covariance matrix for times [0...T-2]
    R : [n_dim_obs, n_dim_obs] array-like
        observation covariance matrix for times [0...T-1]
    b : [T-1, n_dim_state] or [n_dim_state] array-like
        state offsets for times [0...T-2]
    d : [T, n_dim_obs] or [n_dim_obs] array-like
        observation offset for times [0...T-1]
    mu_0 : [n_dim_state] array-like
        mean of initial state distribution
    sigma_0 : [n_dim_state, n_dim_state] array-like
        covariance of initial state distribution
    random_state : optional, numpy random state
        random number generator used in sampling
    em_vars : optional, subset of ['A', 'C', 'b', 'd', 'Q', 'R', 'mu_0',
    'sigma_0'] or 'all'
        if `em_vars` is an iterable of strings only variables in `em_vars`
        will be estimated using EM.  if `em_vars` == 'all', then all
        variables will be estimated.
    """
    def __init__(self, A, C, Q, R, b, d, mu_0, sigma_0, random_state=None,
                 em_vars=['Q', 'R', 'mu_0', 'sigma_0']):
        """Initialize Kalman Filter"""
        self.A = array2d(A)
        self.C = array2d(C)
        self.Q = array2d(Q)
        self.R = array2d(R)
        self.b = array1d(b)
        self.d = array1d(d)
        self.mu_0 = array1d(mu_0)
        self.sigma_0 = array2d(sigma_0)
        self.random_state = random_state
        self.em_vars = em_vars


    def sample(self, T, x_0=None, random_state=None):
        """Sample a trajectory of `T` timesteps in length.

        Parameters
        ==========
        T : int
            number of timesteps

        Returns
        =======
        x : [T, n_dim_state] array
            hidden states corresponding to times [0...T-1]
        z : [T, n_dim_obs] array
            observations corresponding to times [0...T-1]
        """
        n_dim_state = self.A.shape[-2]
        n_dim_obs = self.C.shape[-2]
        x = np.zeros((T, n_dim_state))
        Z = np.zeros((T, n_dim_obs))

        # logic for instantiating rng
        if random_state is None:
            if self.random_state is None:
                rng = np.random.RandomState()
            else:
                rng = self.random_state
        else:
            rng = random_state

        # logic for selecting initial state
        if x_0 is None:
            x_0 = rng.multivariate_normal(self.mu_0, self.sigma_0)

        # logic for generating samples
        for t in range(T):
            if t == 0:
                x[t] = x_0
            else:
                A_t1 = _last_dims(self.A, t - 1)
                b_t1 = _last_dims(self.b, t - 1, ndims=1)
                Q_t1 = _last_dims(self.Q, t - 1)
                x[t] = A_t1.dot(x[t - 1]) + b_t1 \
                    + rng.multivariate_normal(np.zeros(n_dim_state),
                        Q_t1.newbyteorder('='))

            C_t = _last_dims(self.C, t)
            d_t = _last_dims(self.d, t, ndims=1)
            R_t = _last_dims(self.R, t)
            Z[t] = C_t.dot(x[t]) + d_t \
                + rng.multivariate_normal(np.zeros(n_dim_obs),
                    R_t.newbyteorder('='))

        return (x, ma.array(Z))


    def filter(self, Z):
        """Apply the Kalman Filter to estimate the hidden state at time
        t for t = [0...T-1] given observations up to and including time t.
        Observations are assumed to correspond to times [0...T-1].  The output
        of this method corresponding to time T-1 can be used in
        :func:`KalmanFilter.filter_update` for online updating.

        Parameters
        ==========
        Z : [T, n_dim_obs] array-like
            observations corresponding to times [0...T-1].  If `Z` is a masked
            array and any of `Z[t]` is masked, then `Z[t]` will be treated as a
            missing observation.

        Returns
        =======
        mu_filt : [T, n_dim_state]
            mean of hidden state distributions for times [0...T-1] given
            observations up to and including the current time step
        sigma_filt : [T, n_dim_state, n_dim_state] array
            covariance matrix of hidden state distributions for times
            [0...T-1] given observations up to and including the current
            time step
        ll : float
            log likelihood of all observations
        """
        Z = ma.asarray(Z)

        (_, _, _, x_filt, V_filt, ll) = _filter(self.A, self.C,
            self.Q, self.R, self.b, self.d, self.mu_0, self.sigma_0, Z)
        return (x_filt, V_filt, ll)


    def filter_update(self, mu_filt, sigma_filt, z=None, b=None, d=None, t=None):
        """Perform a one-step update to estimate the state at time t+1
        give an observation at time t+1 and the previous estimate for
        time t given observations from times [1...t].  This method is useful if
        one wants to track an object with streaming observations.

        Parameters
        ==========
        mu_filt : [n_dim_state] array
            mean estimate for state at time t given observations from times
            [1...t]
        sigma_filt : [n_dim_state, n_dim_state] array
            covariance of estimate for state at time t given observations from
            times [1...t]
        z : [n_dim_obs] array or None
            observation from time t+1.  If `z` is a masked array and any of
            `z`'s components are masked or if `z` is None, then `z` will be
            treated as a missing observation.
        b : optional, [n_dim_state] array
            state offset for transition from time t to t+1.  If unspecified,
            `self.b` will be used.
        d : optional, [n_dim_obs] array
            observation offset for time t+1.  If unspecified, `self.d` will be
            used.
        t : optional, int
            time step of `mu_filt` and `sigma_filt`.  Used to identify `A`,
            `C`, `b`, and `d` used to transition from time t to t+1 and
            generate observation at time t+1.  If all of the above are constant
            for all time or are specified in this method's arguments, then this
            argument can be disregarded.

        Returns
        =======
        mu_new : [n_dim_state] array
            mean estimate for state at time t+1 given observations from times
            [1...t+1]
        sigma_new : [n_dim_state, n_dim_state] array
            covariance of estimate for state at time t+1 given observations from
            times [1...t+1]
        ll : float
            likelihood of observation z given all previous observations.
        """
        
        if t is not None:
            if b is None:
                b = _last_dims(self.b, t, ndims=1)
            if d is None:
                d = _last_dims(self.d, t + 1, ndims=1)
            A = _last_dims(self.A, t)
            C = _last_dims(self.C, t + 1)
        else:
            if b is None:
                b = self.b
                if len(b.shape) > 1:
                    raise ValueError('b is not constant for all time, you must specify t')
            if d is None:
                d = self.d
                if len(d.shape) > 1:
                    raise ValueError('d is not constant for all time, you must specify t')
            A = self.A
            if len(A.shape) != 2:
                raise ValueError('A is not constant for all time, you must specify t')
            C = self.C
            if len(C.shape) != 2:
                raise ValueError('C is not constant for all time, you must specify t')
        Q = self.Q  # assumed to be time-invariant
        R = self.R

        # Make a masked observation if necessary
        if z is None:
            n_dim_obs = R.shape[0]
            z = ma.array(np.zeros(n_dim_obs))
            z.mask = True

        (mu_pred, sigma_pred) = _filter_predict(A, Q, b, mu_filt, sigma_filt)
        (_, mu_new, sigma_new, ll) = _filter_correct(C, R, d, mu_pred, sigma_pred, z)

        return (mu_new, sigma_new, ll)


    def predict(self, Z):
        """Apply the Kalman Smoother to estimate the hidden state at time
        t for t = [0...T-1] given all observations.  See
        :func:`sklearn.kalman._smooth` for more complex output

        Parameters
        ==========
        Z : [T, n_dim_obs] array-like
            observations corresponding to times [0...T-1].  If `Z` is a masked
            array and any of `Z[t]` is masked, then `Z[t]` will be treated as a
            missing observation.

        Returns
        =======
        mu_smooth : [T, n_dim_state]
            mean of hidden state distributions for times [0...T-1] given
            all observations
        """
        Z = ma.asarray(Z)

        (mu_pred, sigma_pred, _, mu_filt, sigma_filt, ll) = _filter(
            self.A, self.C, self.Q, self.R, self.b, self.d, self.mu_0,
            self.sigma_0, Z)
        (mu_smooth, _, _) = _smooth(
            self.A, mu_filt, sigma_filt, mu_pred, sigma_pred)
        return mu_smooth


    def fit(self, Z, n_iter=10, em_vars=None):
        """Apply the EM algorithm to estimate all parameters specified by
        `em_vars`.  Note that all variables estimated are assumed to be
        constant for all time.  See :func:`sklearn.kalman._em` for details.

        Parameters
        ==========
        Z : [T, n_dim_obs] array-like
            observations corresponding to times [0...T-1].  If `Z` is a masked
            array and any of `Z[t]`'s components is masked, then `Z[t]` will be
            treated as a missing observation.
        n_iter : int, optional
            number of EM iterations to perform
        em_vars : iterable of strings or 'all'
            variables to perform EM over.  Any variable not appearing here is
            left untouched.
        """
        Z = ma.asarray(Z)

        # Create dictionary of variables not to perform EM on
        if em_vars is None:
            em_vars = self.em_vars
        if em_vars == 'all':
            given = {}
        else:
            given = {
                'A': self.A,
                'C': self.C,
                'b': self.b,
                'd': self.d,
                'Q': self.Q,
                'R': self.R,
                'mu_0': self.mu_0,
                'sigma_0': self.sigma_0
            }
            em_vars = set(em_vars)
            for k in given.keys():
                if k in em_vars:
                    given.pop(k)

        # If a parameter is time varying, print a warning
        for (k,v) in self.get_params().items():
            if k in DIM and (not k in given) and len(v.shape) != DIM[k]:
                warn_str = '%s has %s dimensions now; ' % (k, len(v.shape)) + \
                    'after fitting, it will have dimension %d' % (DIM[k],)
                warnings.warn(warn_str)

        # Actual EM iterations
        for i in range(n_iter):
            (mu_pred, sigma_pred, K, mu_filt, sigma_filt, loglik) = _filter(
                self.A, self.C, self.Q, self.R, self.b, self.d, self.mu_0,
                self.sigma_0, Z)
            (mu_smooth, sigma_smooth, L) = _smooth(
                self.A, mu_filt, sigma_filt, mu_pred, sigma_pred)
            sigma_pair_smooth = _smooth_pair(sigma_smooth, L)
            (self.A,  self.C, self.b, self.d, self.Q, self.R, self.mu_0,
                self.sigma_0) = _em( Z, self.b, self.d, mu_smooth,
                    sigma_smooth, sigma_pair_smooth, given=given)
        return self
