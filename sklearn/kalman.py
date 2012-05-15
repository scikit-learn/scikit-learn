'''
Inference for Linear-Gaussian Systems.

This module implements the Kalman Filter, Kalman Smoother, and
EM Algorithm for Linear-Gaussian state space models.    In other
words, assuming that

    x_{t+1} = A_t * x_t + b_t + e_t^1
    y_{t}     = C_t * x_t + d_t + e_t^2
    e_t^1     ~ MultivariateNormal(0, Q)
    e_t^2     ~ MultivariateNormal(0, R)
    x_0         ~ MultivariateNormal(x_0, V_0)

then the Kalman Filter calculates exactly,

    P(x_t | y_{1:t})

the Kalman Smoother,

    P(x_t | y_{1:T})

the EM algorithm,

    argmax_{Q,R} P(y_{1:T}; Q, R)

Observations are assumed to be sampled at times [1...T] while
states are sampled from times [0...T]
'''
import numpy as np


def _last_dims(X, t, ndims=2):
    '''Extract the final `ndim` dimensions at index `t` if
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
    '''
    if len(X.shape) == ndims + 1:
        return X[t]
    elif len(X.shape) == ndims:
        return X
    else:
        raise ValueError("X only has {} ".format(len(X.shape)) +
                         "dimensions when {}".format(ndims) +
                         " or more are required")


def _logmvnpdf(x, mu, sigma):
    '''log density of the multivariate normal distribution

    Parameters
    ==========
    x : [n_dim] array
        point sampled from multivariate normal
    mu : [n_dim] array
        mean of multivariate normal
    sigma : [n_dim, n_dim] array
        covariance of multivariate normal
    '''
    k = len(mu)
    ll = -0.5 * (x - mu).dot(np.linalg.pinv(sigma)).dot(x - mu) + \
        -k * 0.5 * np.log(2 * np.pi) + \
        -0.5 * np.log(np.abs(np.linalg.det(sigma)))
    return ll


def _filter_update(A, C, Q, R, b, d, mu_old, sigma_old, z):
    '''Your usual Kalman update for,

            x'= A*x + b + Normal(0, Q)
            z = C*x + d + Normal(0, R)

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
    mu_new : [n_dim_state] array
        mean of state at time t+1 given observations from times
        [1...t+1]
    sigma_new : [n_dim_state, n_dim_state] array
        covariance of state at time t+1 given observations from times
        [1...t+1]
    ll : float
        log likelihood of observation at time t+1 given observations
        from times [1...t]
    '''
    mu_pred = A.dot(mu_old) + b
    sigma_pred = A.dot(sigma_old).dot(A.T) + Q

    mu_obs = C.dot(mu_pred) + d
    sigma_obs = C.dot(sigma_pred).dot(C.T) + R
    ll = _logmvnpdf(z, mu_obs, sigma_obs)

    K_new = sigma_pred.dot(C.T).dot(np.linalg.pinv(sigma_obs))

    mu_new = mu_pred + K_new.dot(z - mu_obs)
    sigma_new = sigma_pred - K_new.dot(C).dot(sigma_pred)

    return (mu_new, sigma_new, ll)


def _filter(A, C, Q, R, b, d, mu_0, sigma_0, Z):
    '''Apply the Kalman Filter'''
    T = Z.shape[0]
    n_dim_state = len(mu_0)

    x_filt = np.zeros((T + 1, n_dim_state))
    V_filt = np.zeros((T + 1, n_dim_state, n_dim_state))
    ll = np.zeros(T + 1)
    for t in range(T):
        if t == 0:
            x_filt[t] = mu_0
            V_filt[t] = sigma_0
        (x_filt[t + 1], V_filt[t + 1], ll[t + 1]) = _filter_update(
            _last_dims(A, t), _last_dims(C, t), _last_dims(Q, t),
            _last_dims(R, t), _last_dims(b, t, ndims=1),
            _last_dims(d, t, ndims=1), x_filt[t], V_filt[t], Z[t])
    return (x_filt, V_filt, ll)


def _smooth_update(A, Q, b, mu_for, sigma_for, mu_rev, sigma_rev):
    '''One iteration of Kalman Smoothing.    Calculates posterior
    distribution of the hidden state at time t given the observations
    from times [1...T].

    Parameters
    ==========
    A : [n_dim_state, n_dim_state} array
        state transition matrix from time t to t+1
    Q : [n_dim_state, n_dim_state] array
        covariance matrix for state transition from time t to t+1
    b : [n_dim_state] array
        offset for state transition from time t to t+1
    mu_for : [n_dim_state] array
        mean of filtered state at time t given observations from
        times [1...t]
    sigma_for : [n_dim_state, n_dim_state] array
        covariance of filtered state at time t given observations from
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
    '''
    mu_pred = A.dot(mu_for) + b
    sigma_pred = A.dot(sigma_for).dot(A.T) + Q
    L = sigma_for.dot(A.T).dot(np.linalg.pinv(sigma_pred))

    sigma = sigma_for + L.dot(sigma_rev - sigma_pred).dot(L.T)
    mu = mu_for + L.dot(mu_rev - mu_pred)

    return (mu, sigma, L)


def _smooth(A, Q, b, mu_filt, sigma_filt):
    '''Kalman Smoother'''
    T, n_dim_state = mu_filt.shape
    T -= 1    # mu_filt is actually T+1 in length

    mu_smooth = np.zeros((T + 1, n_dim_state))
    sigma_smooth = np.zeros((T + 1, n_dim_state, n_dim_state))
    L = np.zeros((T, n_dim_state, n_dim_state))

    mu_smooth[-1] = mu_filt[-1]
    sigma_smooth[-1] = sigma_filt[-1]

    for t in reversed(range(T)):
        (mu_smooth[t], sigma_smooth[t], L[t]) = _smooth_update(
            _last_dims(A, t), _last_dims(Q, t), _last_dims(b, t, ndims=1),
            mu_filt[t], sigma_filt[t], mu_smooth[t + 1], sigma_smooth[t + 1]
            )
    return (mu_smooth, sigma_smooth, L)


def _em(A, C, b, d, Z, mu_smooth, sigma_smooth, L):
    '''Estimate Q and R given smoothed estimates for states

    Parameters
    ==========
    A : [n_dim_state, n_dim_state} array
        state transition matrix from time t to t+1
    C : [n_dim_obs, n_dim_state] array
        observation matrix for time t+1
    b : [n_dim_state] array
        offset for state transition from time t to t+1
    d : [n_dim_obs] array
        offset for observation at time t+1
    Z : [T, n_dim_obs] array
        observations from times [1...T]
    mu_smooth : [T+1, n_dim_state] array
        means of hidden state for times [0...T] given observations
        from times [1..T]
    sigma_smooth : [T+1, n_dim_state, n_dim_state] array
        covariance of hidden state for times [0...T] given
        observations from times [1...T]
    L : [T, n_dim_state, n_dim_state] array
        gain matrix used in Kalman Smoothing for times [1...T]

    Returns
    =======
    Q : [n_dim_state, n_dim_state] array
        estimated covariance matrix for state transitions
    R : [n_dim_obs, n_dim_obs] array
        estimated covariance matrix for observations
    '''
    T, n_dim_state = mu_smooth.shape
    n_dim_obs = C.shape[-2]
    T -= 1    # mu_smooth is [T+1, n_dim_state] actually

    Q = np.zeros((n_dim_state, n_dim_state))
    for t in range(T):
        A_t = _last_dims(A, t)
        b_t = _last_dims(b, t, ndims=1)
        L_t = L[t]
        x_pred = A_t.dot(mu_smooth[t]) + b_t
        x_next = mu_smooth[t + 1]
        V_now = sigma_smooth[t]
        V_next = sigma_smooth[t + 1]
        Q += np.outer(x_next - x_pred, x_next - x_pred) + A_t.dot(V_now).dot(A_t.T) \
                + V_next - V_next.dot(L_t.T).dot(A_t.T) - A_t.dot(L_t).dot(V_next)
    Q *= 1.0 / T

    R = np.zeros((n_dim_obs, n_dim_obs))
    for t in range(1, T + 1):
        C_t = _last_dims(C, t - 1)
        d_t = _last_dims(d, t - 1, ndims=1)
        y_pred = C_t.dot(mu_smooth[t]) + d_t
        y_t = Z[t - 1]
        V_now = sigma_smooth[t]
        R += np.outer(y_t - y_pred, y_t - y_pred) + C_t.dot(V_now).dot(C_t.T)
    R *= 1.0 / T

    return (Q, R)


class KalmanFilter(object):
    '''
    Implements the Kalman Filter, Kalman Smoother, and EM algorithm for
    a Linear Gaussian model specified by,

        x_{t+1} = A_{t}*x_{t} + b_{t} + Normal(0, Q_{t})
        z_{t}     = C_{t}*x_{t} + d_{t} + Normal(0, R_{t})

    The Kalman Filter is an algorithm designed to estimate P(x_t | z_{1:t}).
    As all state transitions and observations are linear with Gaussian
    distributed noise, these distributions can be represented exactly as
    Gaussian distributions with mean `mu_filt[t]` and covariances
    `sigma_filt[t]`.

    Similarly, the Kalman Smoother is an algorithm designed to estimate
    P(x_t | z_{1:T}).

    The EM algorithm aims to find
        max_{Q,R} P(z_{1:T}; Q,R)
    If we define L(x_{1:t}, Q,R) = log P(z_{1:T}, x_{1:T}; Q, R), then
    the EM algorithm works by iteratively finding,
        P(x_{1:T} | z_{1:T}, Q_{i}, R_{i})
    then by maximizing,
        Q_{i+1}, R_{i+1} = argmax_{Q,R}
            E_{x:1:T} [ L(x_{1:t}, Q,R); z_{1:T], Q_{i}, R_{i} ]
    '''
    def __init__(self, A, C, Q, R, b, d, x_0, V_0,
                             rng=np.random.RandomState(0)):
        '''Initialize Kalman Filter

        Parameters
        ==========
        A : [T, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
            state transition matrix
        C : [T, n_dim_obs, n_dim_obs] or [n_dim_obs, n_dim_obs] array-like
            observation matrix
        Q : [T, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
            state transition covariance matrix
        R : [T, n_dim_state, n_dim_state] or [n_dim_state, n_dim_state] array-like
            observation covariance matrix
        b : [T, n_dim_state] or [n_dim_state] array-like
            state offset
        d : [T, n_dim_obs] or [n_dim_obs] array-like
            observation offset
        x_0 : [n_dim_state] array-like
            mean of initial state distribution
        V_0 : [n_dim_state, n_dim_state] array-like
            covariance of initial state distribution
        rng : numpy random state
            random number generator
        '''
        self.A = np.asarray(A)
        self.C = np.asarray(C)
        self.Q = np.asarray(Q)
        self.R = np.asarray(R)
        self.b = np.asarray(b)
        self.d = np.asarray(d)
        self.x_0 = np.asarray(x_0)
        self.V_0 = np.asarray(V_0)
        self.rng = rng

    def sample(self, T, x_0=None):
        '''Sample a trajectory of T timesteps in length.

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
        '''
        n_dim_state = self.A.shape[-2]
        n_dim_obs = self.C.shape[-2]
        x = np.zeros((T + 1, n_dim_state))
        z = np.zeros((T, n_dim_obs))

        if x_0 is None:
            x_0 = self.x_0

        for t in range(T + 1):
            if t == 0:
                x[t] = x_0
            else:
                x[t] = self.A.dot(x[t - 1]) + _last_dims(self.b, t, ndims=1) \
                    + self.rng.multivariate_normal(np.zeros(n_dim_state),
                                 _last_dims(self.Q, t).newbyteorder('='))
                z[t - 1] = self.C.dot(x[t]) + _last_dims(self.d, t, ndims=1) \
                    + self.rng.multivariate_normal(np.zeros(n_dim_obs),
                        _last_dims(self.R, t).newbyteorder('='))
        return (x, z)

    def filter(self, Z):
        '''Apply the Kalman Filter to estimate the hidden state at time
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
        '''
        Z = np.asarray(Z)
        return _filter(self.A, self.C, self.Q, self.R, self.b, self.d,
                       self.x_0, self.V_0, Z)

    def smooth(self, Z):
        '''Apply the Kalman Smoother to estimate the hidden state at time
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
        '''
        Z = np.asarray(Z)
        (mu_filt, sigma_filt, ll) = self.filter(Z)
        (mu_smooth, sigma_smooth, _) = _smooth(self.A, self.Q, self.b,
                                               mu_filt, sigma_filt)
        return (mu_smooth, sigma_smooth, ll)

    def em(self, Z, n_iter=10):
        '''Apply the EM algorithm to estimate the covariance matrices
        `Q` and `R` representing the covariance for state transitions
        and observations, respectively.

        Parameters
        ==========
        Z : [T, n_dim_obs] array-like
            observations corresponding to times [1...T]
        n_iter : int, optional
            number of EM iterations to perform

        Returns
        =======
        Q : [n_dim_state, n_dim_state] array
            estimated covariance of state transitions
        R : [n_dim_obs, n_dim_obs] array
            estimated covariance of observations
        ll : float
            log likelihood of all observations
        '''
        Z = np.asarray(Z)
        ll = np.zeros(n_iter)
        for i in range(n_iter):
            (mu_filt, sigma_filt, loglik) = self.filter(Z)
            (mu_smooth, sigma_smooth, L) = _smooth(self.A, self.Q, self.b,
                                                   mu_filt, sigma_filt)
            (self.Q, self.R) = _em(self.A, self.C, self.b, self.d, Z, mu_smooth,
                                   sigma_smooth, L)
            ll[i] = np.sum(loglik)
        return (self.Q, self.R, ll)


if __name__ == '__main__':
    from pprint import pprint
    import matplotlib.pyplot as plt
    from sklearn.datasets import load_kalman_data
    data = load_kalman_data()

    kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q, R=data.R, b=data.b,
                                        d=data.d, x_0=data.x_0, V_0=data.V_0)
    (x, z) = kf.sample(100)
    assert x.shape == (101, data.A.shape[0])
    assert z.shape == (100, data.C.shape[0])
    print 'Sampling OK'

    (x_filt, V_filt, _) = kf.filter(Z=data.data)
    print 'Filtering OK'

    (x_smooth, V_smooth, _) = kf.smooth(Z=data.data)
    print 'Smoothing OK'

    kf.Q = 10 * np.eye(5)
    kf.R = 10 * np.eye(2)
    (Q, R, ll) = kf.em(Z=data.data, n_iter=10)
    print 'EM OK'

    plt.figure()
    plt.hold(True)
    lines_true = plt.plot(data.target, linestyle='-', label='true state')
    lines_filt = plt.plot(x_filt, linestyle='--', color='g', label='filtered state')
    lines_smooth = plt.plot(x_smooth, linestyle='-.', color='r', label='smoothed state')
    plt.legend((lines_true[0], lines_filt[0], lines_smooth[0]),
                ('true state', 'filtered state', 'smoothed state'))
    plt.xlabel('time')
    plt.ylabel('state')

    plt.figure()
    plt.plot(ll)
    plt.xlabel('iteration number')
    plt.ylabel('log likelihood')
    plt.show()
