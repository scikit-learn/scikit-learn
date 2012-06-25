'''
==============================================
Online State Estimation with the Kalman Filter
==============================================

The Kalman Filter updates the state mean and covariance matrix in a recursive
fashion and is thus ideal for online state estimation.  This example shows how
it can be applied with the :mod:`sklearn.kalman` module.
'''
import numpy as np
import pylab as pl

from sklearn.datasets import load_kalman_data
from sklearn.kalman import KalmanFilter

# Initialize the Kalman Filter
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
    random_state=0
)

# Estimate mean and covariance of hidden state distribution iteratively.
T = data.data.shape[0]
n_dim_state = data.transition_matrix.shape[0]
filtered_state_means = np.zeros((T, n_dim_state))
filtered_state_covariances = np.zeros((T, n_dim_state, n_dim_state))
for t in range(T - 1):
    if t == 0:
        filtered_state_means[t] = data.initial_state_mean
        filtered_state_covariances[t] = data.initial_state_covariance
    (filtered_state_means[t + 1],
     filtered_state_covariances[t + 1], _) = kf.filter_update(
         filtered_state_means[t], filtered_state_covariances[t],
         data.data[t + 1], transition_offset=data.transition_offsets[t],
         observation_offset=data.observation_offset
    )

# draw estimates
pl.figure()
pl.hold(True)
lines_true = pl.plot(data.target, color='b')
lines_filt = pl.plot(filtered_state_means, color='r')
pl.legend((lines_true[0], lines_filt[0]), ('true', 'filt'))
pl.show()
