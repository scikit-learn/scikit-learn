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

    P(x_t | z_{1:t})

the Kalman Smoother,

.. math::

    P(x_t | z_{1:T})

the EM algorithm, for :math:`\theta = (A,C,Q,R,\mu_0,\Sigma_0)`,

.. math::

    \arg\max_{\theta} P(z_{1:T}; \theta)

Observations are assumed to be sampled at times [1...T] while
states are sampled from times [0...T]

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
import numpy as np

from sklearn.utils import array2d


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
    if len(X.shape) == ndims + 1:
        return X[t]
    elif len(X.shape) == ndims:
        return X
    else:
        raise ValueError("X only has {} ".format(len(X.shape)) +
                         "dimensions when {}".format(ndims) +
                         " or more are required")


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
    k = len(mu)
    ll = -0.5 * (x - mu).dot(np.linalg.pinv(sigma)).dot(x - mu) + \
        -k * 0.5 * np.log(2 * np.pi) + \
        -0.5 * np.log(np.abs(np.linalg.det(sigma)))
    return ll


def _filter_update(A, C, Q, R, b, d, mu_old, sigma_old, z):
    """Your usual Kalman update for,

    .. math::

        x'  &= Ax + b + Normal(0, Q) \\
        z   &= Cx + d + Normal(0, R)

    Calculates the posterior distribution over the hidden state at
    time t+1 given observations from times [1...t+1]

    Parameters
    ==========
    A : [n_dim_state, n_dim_state} array
        state transition matrix from time t to t+1
    C : [n_dim_obs, n_dim_state] array
        observation matrix for time t+1
    Q : [n_dim_state, n_dim_state] array
        covariance matrix for state transition from time t to t+1
    R : [n_dim_obs, n_dim_obs] array
        covariance matrix for observation at time t+1
    b : [n_dim_state] array
        offset for state transition from time t to t+1
    d : [n_dim_obs] array
        offset for observation at time t+1
    mu_old : [n_dim_state]
        mean of state at time t given observations from times
        [1...t]
    sigma_old : [n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [1...t]
    z : [n_dim_obs] array
        observation at time t+1

    Returns
    =======
    mu_pred : [n_dim_state] array
        mean of state at time t+1 given observations from times [1...t]
    sigma_pred : [n_dim_state, n_dim_state] array
        covariance of state at time t+1 given observations from times
        [1...t]
    K : [n_dim_state, n_dim_obs] array
        Kalman gain matrix for time t+1
    mu_new : [n_dim_state] array
        mean of state at time t+1 given observations from times
        [1...t+1]
    sigma_new : [n_dim_state, n_dim_state] array
        covariance of state at time t+1 given observations from times
        [1...t+1]
    ll : float
        log likelihood of observation at time t+1 given observations
        from times [1...t]
    """
    mu_pred = A.dot(mu_old) + b
    sigma_pred = A.dot(sigma_old).dot(A.T) + Q

    mu_obs = C.dot(mu_pred) + d
    sigma_obs = C.dot(sigma_pred).dot(C.T) + R
    ll = _logmvnpdf(z, mu_obs, sigma_obs)

    K_new = sigma_pred.dot(C.T).dot(np.linalg.pinv(sigma_obs))

    mu_new = mu_pred + K_new.dot(z - mu_obs)
    sigma_new = sigma_pred - K_new.dot(C).dot(sigma_pred)

    return (mu_pred, sigma_pred, K_new, mu_new, sigma_new, ll)


def _filter(A, C, Q, R, b, d, mu_0, sigma_0, Z):
    """Apply the Kalman Filter.  See :func:`KalmanFilter.filter` for
    details.
    """
    T = Z.shape[0] - 1
    n_dim_state = len(mu_0)
    n_dim_obs = Z.shape[1]

    x_pred = np.zeros((T + 1, n_dim_state))
    V_pred = np.zeros((T + 1, n_dim_state, n_dim_state))
    K = np.zeros((T + 1, n_dim_state, n_dim_obs))
    x_filt = np.zeros((T + 1, n_dim_state))
    V_filt = np.zeros((T + 1, n_dim_state, n_dim_state))
    ll = np.zeros(T + 1)
    for t in range(T):
        if t == 0:
            x_filt[t] = x_pred[t] = mu_0
            V_filt[t] = V_pred[t] = sigma_0
        (x_pred[t + 1], V_pred[t + 1], K[t + 1],
            x_filt[t + 1], V_filt[t + 1], ll[t + 1]) = _filter_update(
                _last_dims(A, t), _last_dims(C, t + 1), _last_dims(Q, t),
                _last_dims(R, t + 1), _last_dims(b, t, ndims=1),
                _last_dims(d, t + 1, ndims=1), x_filt[t], V_filt[t], Z[t + 1])
    return (x_pred, V_pred, K, x_filt, V_filt, ll)


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
    """Apply the Kalman Smoother.  See :func:`KalmanFilter.smooth` for
    details.
    """
    T, n_dim_state = mu_filt.shape
    T -= 1    # mu_filt is actually T+1 in length

    mu_smooth = np.zeros((T + 1, n_dim_state))
    sigma_smooth = np.zeros((T + 1, n_dim_state, n_dim_state))
    L = np.zeros((T + 1, n_dim_state, n_dim_state))

    mu_smooth[-1] = mu_filt[-1]
    sigma_smooth[-1] = sigma_filt[-1]

    for t in reversed(range(T)):
        (mu_smooth[t], sigma_smooth[t], L[t]) = _smooth_update(
            _last_dims(A, t),
            mu_filt[t], sigma_filt[t],
            mu_pred[t + 1], sigma_pred[t + 1],
            mu_smooth[t + 1], sigma_smooth[t + 1],
            )
    return (mu_smooth, sigma_smooth, L)


def _smooth_pair(sigma_smooth, L):
    """
    sigma_pair_smooth[t] = :math:`Cov(x_t, x_{t-1} | Y_{1:T})`
    """
    T, n_dim_state, _ = sigma_smooth.shape
    T -= 1
    sigma_pair_smooth = np.zeros((T + 1, n_dim_state, n_dim_state))
    for t in range(1, T + 1):
        sigma_pair_smooth[t] = sigma_smooth[t].dot(L[t - 1].T)
    return sigma_pair_smooth


def _em(Z, b, d, mu_smooth, sigma_smooth, sigma_pair_smooth, L, given={}):
    """Estimate Linear-Gaussian model parameters.

    Parameters
    ==========
    Z : [T, n_dim_obs] array
        observations for times [0...T] (t=0 is ignored)
    b : [n_dim_state] or [T, n_dim_state] array
        transition offset
    d : [n_dim_obs] or [T, n_dim_obs] array
        observation offsets
    mu_smooth : [T+1, n_dim_state] array
        mu_smooth[t] = mean of state at time t given all observations
    sigma_smooth : [T+1, n_dim_state, n_dim_state] array
        sigma_smooth[t] = covariance of state at time t given all observations
    sigma_pair_smooth : [T+1, n_dim_state, n_dim_state] array
        sigma_pair_smooth[t] = covariance between states at times t and
        t-1 given all observations
    given: dict
        if one of the variables EM is capable of predicting is in given, then
        that value will be used and EM will not attempt to estimate it.  e.g.,
        if 'C' is in given, C will not be estimated and give['C'] will be
        returned in its place.

    Returns
    =======
    A : [n_dim_state, n_dim_state] array
        estimated transition matrix
    C : [n_dim_obs, n_dim_state] array
        estimated observation matrix
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

    return (A, C, Q, R, mu_0, sigma_0)


def _em_C(Z, d, mu_smooth, sigma_smooth):
    """Maximize expected log likelihood of observations with respect to the
    observation matrix C.

    .. math::

        C &= ( \sum_{t=1}^{T} (y_t - d_t) \E[x_t] )
             ( \sum_{t=1}^{T} E[x_t x_t^T] )^-1
    """
    _, n_dim_state = mu_smooth.shape
    T, n_dim_obs = Z.shape
    T -= 1
    res1 = np.zeros((n_dim_obs, n_dim_state))
    res2 = np.zeros((n_dim_state, n_dim_state))
    for t in range(1, T + 1):
        d_t = _last_dims(d, t, ndims=1)
        res1 += np.outer(Z[t] - d_t, mu_smooth[t])
        res2 += sigma_smooth[t] + np.outer(mu_smooth[t], mu_smooth[t])
    return res1.dot(np.linalg.pinv(res2))


def _em_R(Z, d, C, mu_smooth, sigma_smooth):
    """Maximize expected log likelihood of observations with respect to the
    observation covariance matrix R.

    .. math::

        R &= \frac{1}{T} \sum_{t=1}^T
               [y_t - C \E[x_t] - d_t] [y_t - C \E[x_t] - d_t]^T
               + C Var(x_t) C^T
    """
    _, n_dim_state = mu_smooth.shape
    T, n_dim_obs = Z.shape
    T -= 1
    res = np.zeros((n_dim_obs, n_dim_obs))
    for t in range(1, T + 1):
        C_t = _last_dims(C, t)
        d_t = _last_dims(d, t, ndims=1)
        err = Z[t] - C_t.dot(mu_smooth[t]) - d_t
        res += np.outer(err, err) \
            + C_t.dot(sigma_smooth[t]).dot(C_t.T)
    return (1.0 / T) * res


def _em_A(b, mu_smooth, sigma_smooth, sigma_pair_smooth):
    """Maximize expected log likelihood of observations with respect to the
    transition matrix A.

    .. math::

        A &= ( \sum_{t=1}^{T} \E[x_t x_{t-1}^{T}] - b_{t-1} \E[x_{t-1}]^T )
             ( \sum_{t=1}^{T} \E[x_{t-1} x_{t-1}^T] )^{-1}
    """
    T, n_dim_state, _ = sigma_smooth.shape
    T -= 1
    res1 = np.zeros((n_dim_state, n_dim_state))
    res2 = np.zeros((n_dim_state, n_dim_state))
    for t in range(1, T + 1):
        res1 += sigma_smooth[t - 1] \
            + np.outer(mu_smooth[t - 1], mu_smooth[t - 1])
        res2 += sigma_pair_smooth[t] \
            + np.outer(mu_smooth[t], mu_smooth[t - 1]) \
            - np.outer(_last_dims(b, t - 1, ndims=1), mu_smooth[t - 1])
    return res2.dot(np.linalg.pinv(res1))


def _em_Q(A, b, mu_smooth, sigma_smooth, sigma_pair_smooth):
    """Maximize expected log likelihood of observations with respect to the
    transition covariance matrix Q.

    .. math::

        Q &= \frac{1}{T} \sum_{t=0}^{T-1}
               (\E[x_{t+1}] A \E[x_t] - b_t) (\E[x_{t+1}] A \E[x_t] - b_t)^T
               + A Var(x_t) A^T + Var(x_{t+1}) - Var(x_{t+1}) L_t^T A^T
               - A_t L_t Var(x_{t+1})
    """
    T, n_dim_state, _ = sigma_smooth.shape
    T -= 1
    res = np.zeros((n_dim_state, n_dim_state))
    for t in range(T):
        A_t = _last_dims(A, t)
        b_t = _last_dims(b, t, ndims=1)
        err = mu_smooth[t + 1] - (A_t.dot(mu_smooth[t]) + b_t)
        Vt1t_A = sigma_pair_smooth[t + 1].dot(A_t.T)
        res += np.outer(err, err) + A_t.dot(sigma_smooth[t]).dot(A_t.T) \
            + sigma_smooth[t + 1] - Vt1t_A - Vt1t_A.T

    return (1.0 / T) * res


def _em_mu_0(mu_smooth):
    """Maximize expected log likelihood of observations with respect to the
    initial state distribution mean mu_0.
    """

    return mu_smooth[0]


def _em_sigma_0(mu_0, mu_smooth, sigma_smooth):
    """Maximize expected log likelihood of observations with respect to the
    covariance of the initial state distribution sigma_0
    """
    x0 = mu_smooth[0]
    x0_x0 = sigma_smooth[0] + np.outer(x0, x0)
    return x0_x0 - np.outer(mu_0, x0) - np.outer(x0, mu_0) \
        + np.outer(mu_0, mu_0)


class KalmanFilter(object):
    """
    Implements the Kalman Filter, Kalman Smoother, and EM algorithm for
    a Linear Gaussian model specified by,

    .. math::

        x_{t+1}   &= A_{t} x_{t} + b_{t} + Normal(0, Q_{t}) \\
        z_{t}     &= C_{t} x_{t} + d_{t} + Normal(0, R_{t})

    The Kalman Filter is an algorithm designed to estimate 
    :math:`P(x_t | z_{1:t})``.  As all state transitions and observations are
    linear with Gaussian distributed noise, these distributions can be
    represented exactly as Gaussian distributions with mean `mu_filt[t]` and
    covariances `sigma_filt[t]`.

    Similarly, the Kalman Smoother is an algorithm designed to estimate
    :math:`P(x_t | z_{1:T})`.

    The EM algorithm aims to find for :math:`theta = (A,C,Q,R,\mu_0,\sigma_0)`

    .. math::

        \max_{\theta} P(z_{1:T}; \theta)

    If we define :math:`L(x_{1:t},\theta) = \log P(z_{1:T}, x_{1:T}; \theta)`,
    then the EM algorithm works by iteratively finding,

    .. math::

        P(x_{1:T} | z_{1:T}, \theta_i)

    then by maximizing,

    .. math::

        \theta_{i+1} = \arg\max_{\theta}
            E_{x_{1:T}} [ L(x_{1:t}, \theta)| z_{1:T], \theta_i ]

    Parameters
    ==========
    A : [T, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
        state transition matrix
    C : [T, n_dim_obs, n_dim_obs] or [n_dim_obs, n_dim_obs] array-like
        observation matrix
    Q : [T, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
        state transition covariance matrix
    R : [T, n_dim_obs, n_dim_obs] or [n_dim_obs, n_dim_obs] array-like
        observation covariance matrix
    b : [T, n_dim_state] or [n_dim_state] array-like
        state offset
    d : [T, n_dim_obs] or [n_dim_obs] array-like
        observation offset
    x_0 : [n_dim_state] array-like
        mean of initial state distribution
    V_0 : [n_dim_state, n_dim_state] array-like
        covariance of initial state distribution
    rng : optional, numpy random state
        random number generator used in sampling
    em_vars : optional, subset of ['A', 'C', 'Q', 'R', 'x_0', 'V_0'] or 'all'
        if `em_vars` is an iterable of strings only variables in `em_vars`
        will be estimated using EM.  if `em_vars` == 'all', then all
        variables will be estimated.
    """
    def __init__(self, A, C, Q, R, b, d, x_0, V_0, random_state=None,
                 em_vars='all'):
        """Initialize Kalman Filter"""
        self.A = array2d(A)
        self.C = array2d(C)
        self.Q = array2d(Q)
        self.R = array2d(R)
        self.b = np.asarray(np.atleast_1d(b))
        self.d = np.asarray(np.atleast_1d(d))
        self.x_0 = np.asarray(np.atleast_1d(x_0))
        self.V_0 = array2d(V_0)
        self.random_state = random_state
        self.em_vars = em_vars

        # TODO make this automatic
        if len(self.R.shape) > 2:
            self.R = _pad(self.R)
        if len(self.C.shape) > 2:
            self.C = _pad(self.C)
        if len(self.d.shape) > 1:
            self.d = _pad(self.d)


    def sample(self, T, x_0=None):
        """Sample a trajectory of `T` timesteps in length.

        Parameters
        ==========
        T : int
            number of timesteps

        Returns
        =======
        x : [T+1, n_dim_state] array
            hidden states corresponding to times [0...T]
        z : [T, n_dim_obs] array
            observations corresponding to times [1...T]
        """
        n_dim_state = self.A.shape[-2]
        n_dim_obs = self.C.shape[-2]
        x = np.zeros((T + 1, n_dim_state))
        Z = np.zeros((T + 1, n_dim_obs))

        if x_0 is None:
            x_0 = self.x_0

        if self.random_state is None:
            rng = np.random.RandomState()
        else:
            rng = self.random_state

        for t in range(T + 1):
            if t == 0:
                x[t] = x_0
            else:
                A_t = _last_dims(self.A, t)
                b_t = _last_dims(self.b, t, ndims=1)
                C_t = _last_dims(self.C, t)
                d_t = _last_dims(self.d, t, ndims=1)
                Q_t = _last_dims(self.Q, t)
                R_t = _last_dims(self.R, t)
                x[t] = A_t.dot(x[t - 1]) + b_t \
                    + rng.multivariate_normal(np.zeros(n_dim_state),
                        Q_t.newbyteorder('='))
                Z[t] = C_t.dot(x[t]) + d_t \
                    + rng.multivariate_normal(np.zeros(n_dim_obs),
                        R_t.newbyteorder('='))
        return (x, Z[1:])


    def filter(self, Z):
        """Apply the Kalman Filter to estimate the hidden state at time
        t for t = [0...T] given observations up to and including time t.
        Observations are assumed to correspond to times [1...T].

        Parameters
        ==========
        Z : [T, n_dim_obs] array-like
            observations corresponding to times [1...T]

        Returns
        =======
        mu_filt : [T+1, n_dim_state]
            mean of hidden state distributions for times [0...T] given
            observations up to and including the current time step
        sigma_filt : [T+1, n_dim_state, n_dim_state] array
            covariance matrix of hidden state distributions for times
            [0...T] given observations up to and including the current
            time step
        ll : float
            log likelihood of all observations
        """
        # Insert empty observation for t=0
        Z = _pad(np.asarray(Z))

        (_, _, _, x_filt, V_filt, ll) = _filter(self.A, self.C,
            self.Q, self.R, self.b, self.d, self.x_0, self.V_0, Z)
        return (x_filt, V_filt, ll)


    def filter_update(self, mu_filt, sigma_filt, z, b=None, d=None, t=None):
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
        z : [n_dim_obs] array
            observation from time t+1
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
                if len(self.b.shape) > 1:
                    raise ValueError('b is not constant for all time, you must specify t')
                else:
                    b = self.b
            if d is None:
                if len(self.d.shape) > 1:
                    raise ValueError('d is not constant for all time, you must specify t')
                else:
                    d = self.d
            A = self.A
            C = self.C
            if len(A.shape) != 2:
                raise ValueError('A is not constant for all time, you must specify t')
            if len(C.shape) != 2:
                raise ValueError('C is not constant for all time, you must specify t')
        Q = self.Q
        R = self.R

        (_, _, _, mu_new, sigma_new, ll) = _filter_update(
            A, C, Q, R, b, d, mu_filt, sigma_filt, z
        )

        return (mu_new, sigma_new, ll)


    def smooth(self, Z):
        """Apply the Kalman Smoother to estimate the hidden state at time
        t for t = [0...T] given all observations.

        Parameters
        ==========
        Z : [T, n_dim_obs] array-like
            observations corresponding to times [1...T]

        Returns
        =======
        mu_smooth : [T+1, n_dim_state]
            mean of hidden state distributions for times [0...T] given
            all observations
        sigma_smooth : [T+1, n_dim_state, n_dim_state] array
            covariance matrix of hidden state distributions for times
            [0...T] given all observations
        ll : float
            log likelihood of all observations
        """
        # Insert empty observation for t=0
        Z = _pad(np.asarray(Z))

        (mu_pred, sigma_pred, _, mu_filt, sigma_filt, ll) = _filter(
            self.A, self.C, self.Q, self.R, self.b, self.d, self.x_0, self.V_0, Z)
        (mu_smooth, sigma_smooth, _) = _smooth(
            self.A, mu_filt, sigma_filt, mu_pred, sigma_pred)
        return (mu_smooth, sigma_smooth, ll)


    def em(self, Z, n_iter=10, em_vars=None):
        """Apply the EM algorithm to estimate all parameters specified by
        `em_vars`.  Note that all variables estimated are assumed to be
        constant for all time.

        Parameters
        ==========
        Z : [T, n_dim_obs] array-like
            observations corresponding to times [1...T]
        n_iter : int, optional
            number of EM iterations to perform
        em_vars : iterable of strings or 'all'
            variables to perform EM over.  Any variable not appearing here is
            left untouched.

        Returns
        =======
        A : [n_dim_state, n_dim_state] array
            estimated transition matrix
        C : [n_dim_obs, n_dim_state] array
            estimated observation matrix
        Q : [n_dim_state, n_dim_state] array
            estimated covariance matrix for state transitions
        R : [n_dim_obs, n_dim_obs] array
            estimated covariance matrix for observations
        mu_0 : [n_dim_state] array
            estimated mean of initial state distribution
        sigma_0 : [n_dim_state] array
            estimated covariance of initial state distribution
        ll : float
            log likelihood of all observations
        """
        # Insert empty observation for t=0
        Z = _pad(np.asarray(Z))

        # Create dictionary of variables not to perform EM on
        if em_vars is None:
            em_vars = self.em_vars
        if em_vars == 'all':
            given = {}
        else:
            given = {
                'A': self.A,
                'C': self.C,
                'Q': self.Q,
                'R': self.R,
                'mu_0': self.x_0,
                'sigma_0': self.V_0
            }
            em_vars = set(em_vars)
            for k in given.keys():
                if k in em_vars:
                    given.pop(k)

        ll = np.zeros(n_iter)
        for i in range(n_iter):
            (mu_pred, sigma_pred, K, mu_filt, sigma_filt, loglik) = _filter(
                self.A, self.C, self.Q, self.R, self.b, self.d, self.x_0, self.V_0, Z)
            (mu_smooth, sigma_smooth, L) = _smooth(
                self.A, mu_filt, sigma_filt, mu_pred, sigma_pred)
            sigma_pair_smooth = _smooth_pair(sigma_smooth, L)
            (self.A,  self.C, self.Q, self.R, self.x_0, self.V_0) = _em(
                Z, self.b, self.d, mu_smooth, sigma_smooth, sigma_pair_smooth, L, given=given)
            ll[i] = np.sum(loglik)
            #print 'likelihood @ iter={}: {}'.format(i, ll[i])
        return (self.A, self.C, self.Q, self.R, self.x_0, self.V_0, ll)


    def fit(self, Z, n_iter=10, em_vars=None):
        """Apply the EM algorithm.  See :func:`KalmanFilter.em` for details."""
        self.em(Z, n_iter=n_iter, em_vars=em_vars)
        return self


    def predict(self, Z):
        """Apply the Kalman Smoother.  See :func:`KalmanFilter.smooth` for details."""
        return self.smooth(Z)
