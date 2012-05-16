'''
===========
Plot Kalman
===========

This example shows how one can use the Kalman Filter, Kalman Smoother,
and EM algorithm to estimate the hidden state of a Linear-Gaussian
dynamics model.  In this model, one assumes that the hidden state
of a discrete-time system is distributed according to a Multivariate
Gaussian distribution at time 0, then evolves like so,

  x_{t+1} = A_t * x_t + b_t + e_t^1
  y_{t}   = C_t * x_t + d_t + e_t^2
  e_t^1   ~ MultivariateNormal(0, Q)
  e_t^2   ~ MultivariateNormal(0, R)
  x_0     ~ MultivariateNormal(x_0, V_0)

The Kalman Filter and Smoother are two methods for estimating the
hidden state x_t given observations y_t.  The EM algorithm, in
addition, allows one to estimate Q and R in an iterative fashion.
'''
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import load_kalman_data
from sklearn.kalman import KalmanFilter

# Load data and initialize Kalman Filter
data = load_kalman_data()

kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q_0, R=data.R_0, b=data.b,
                  d=data.d, x_0=data.x_0, V_0=data.V_0)

# Learn good values for the state transition covariance matrix Q and
# observation covariance matrix R using the EM algorithm.
(Q, R, ll) = kf.em(Z=data.data, n_iter=10)
kf.Q = Q
kf.R = R

# Estimate the state without using an observations.  This will let us see how
# good we could do if we ran blind.
n_dim_state = data.A.shape[0]
T = data.data.shape[0]
x_blind = np.zeros( (T+1, n_dim_state) )
for t in range(T):
  if t == 0:
    x_blind[t] = data.x_0
  x_blind[t+1] = data.A.dot(x_blind[t]) + data.b[t]

# Estimate the hidden states using observations up to and including
# time t for t in [0...T].  This method outputs the mean and covariance
# characterizing the Multivariate Normal distribution for
#   P(x_t | z_{1:t})
(x_filt, V_filt, _) = kf.filter(data.data)

# Estimate the hidden states using all observations.  These estimates
# will be 'smoother' (and are to be preferred) to those produced by
# simply filtering as they are made with later observations in mind.
# Probabilistically, this method produces the mean and covariance
# characterizing,
#    P(x_t | z_{1:T})
(x_smooth, V_smooth, _) = kf.smooth(data.data)

# Draw the true, filtered, and smoothed state estimates for all
# 5 dimensions.
plt.figure()
plt.hold(True)
lines_true = plt.plot(data.target, linestyle='-', color='b')
lines_blind = plt.plot(x_blind, linestyle=':', color='m')
lines_filt = plt.plot(x_filt, linestyle='--', color='g')
lines_smooth = plt.plot(x_smooth, linestyle='-.', color='r')
plt.legend((lines_true[0], lines_blind[0], lines_filt[0], lines_smooth[0]),
            ('true', 'blind', 'filtered', 'smoothed'))
plt.xlabel('time')
plt.ylabel('state')

# Draw the log likelihood of the data as produced by the EM algorithm.
# Notice that it is increasing in the number of iterations.
plt.figure()
plt.plot(ll)
plt.xlabel('iteration number')
plt.ylabel('log likelihood')
plt.show()
