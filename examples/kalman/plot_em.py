'''
=============================
EM for Linear-Gaussian Models
=============================

This example shows how one may use the EM algorithm to estimate model
parameters in a Linear-Gaussian model.
'''
import numpy as np
import pylab as pl

from sklearn.datasets import load_kalman_data
from sklearn.kalman import KalmanFilter

# Load data and initialize Kalman Filter
data = load_kalman_data()
kf = KalmanFilter(
    data.transition_matrix,
    data.observation_matrix,
    data.initial_transition_covariance,
    data.initial_observation_covariance,
    data.transition_offsets,
    data.observation_offset,
    data.initial_state_mean,
    data.initial_state_covariance,
    em_vars=[
      'transition_matrices', 'observation_matrices',
      'transition_covariance', 'observation_covariance',
      'observation_offsets', 'initial_state_mean',
      'initial_state_covariance'
    ])

# Learn good values for variables named in `em_vars` using the EM algorithm
loglikelihoods = np.zeros(10)
for i in range(len(loglikelihoods)):
    kf = kf.fit(X=data.data, n_iter=1)
    loglikelihoods[i] = np.sum(kf.filter(X=data.data)[-1])

# Estimate the state without using any observations.  This will let us see how
# good we could do if we ran blind.
n_dim_state = data.transition_matrix.shape[0]
T = data.data.shape[0]
blind_state_estimates = np.zeros((T, n_dim_state))
for t in range(T - 1):
    if t == 0:
        blind_state_estimates[t] = data.initial_state_mean
    blind_state_estimates[t + 1] = data.transition_matrix.  \
        dot(blind_state_estimates[t]) +   \
        data.transition_offsets[t]

# Estimate the hidden states using observations up to and including
# time t for t in [0...T-1].  This method outputs the mean and covariance
# characterizing the Multivariate Normal distribution for
#   P(x_t | z_{1:t})
(filtered_state_estimates, _, _) = kf.filter(data.data)

# Estimate the hidden states using all observations.  These estimates
# will be 'smoother' (and are to be preferred) to those produced by
# simply filtering as they are made with later observations in mind.
# Probabilistically, this method produces the mean and covariance
# characterizing,
#    P(x_t | z_{1:T})
smoothed_state_estimates = kf.predict(data.data)

# Draw the true, filtered, and smoothed state estimates for all
# 5 dimensions.
pl.figure(figsize=(16, 6))
pl.hold(True)
lines_true = pl.plot(data.target, linestyle='-', color='b')
lines_blind = pl.plot(blind_state_estimates, linestyle=':', color='m')
lines_filt = pl.plot(filtered_state_estimates, linestyle='--', color='g')
lines_smooth = pl.plot(smoothed_state_estimates, linestyle='-.', color='r')
pl.legend(
    (lines_true[0], lines_blind[0], lines_filt[0], lines_smooth[0]),
    ('true', 'blind', 'filtered', 'smoothed')
)
pl.xlabel('time')
pl.ylabel('state')
pl.xlim(xmax=500)

# Draw log likelihood of observations as a function of EM iteration number.
# Notice how it is increasing (this is guaranteed by the EM algorithm)
pl.figure()
pl.plot(loglikelihoods)
pl.xlabel('em iteration number')
pl.ylabel('log likelihood')
pl.show()
