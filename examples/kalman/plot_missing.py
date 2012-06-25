'''
====================================================
Applying the Kalman Filter with Missing Observations
====================================================

This example shows how one may apply all of :mod:`sklearn.kalman`'s Kalman
Smoother, even with missing observations.
'''
import numpy as np
import numpy.ma as ma
import pylab as pl
from sklearn.kalman import KalmanFilter

# specify parameters
random_state = np.random.RandomState(0)
transition_matrix = [[1, 0.1], [0, 1]]
transition_offset = [-0.1, 0.1]
observation_matrix = np.eye(2) + random_state.randn(2, 2) * 0.1
observation_offset = [1.0, -1.0]
transition_covariance = np.eye(2)
observation_covariance = np.eye(2) + random_state.randn(2, 2) * 0.1
initial_state_mean = [5, -5]
initial_state_covariance = [[1, 0.1], [-0.1, 1]]
T = 50

# sample from model
kf = KalmanFilter(
    transition_matrix, observation_matrix, transition_covariance,
    observation_covariance, transition_offset, observation_offset,
    initial_state_mean, initial_state_covariance, random_state=0
)
(states, observations_all) = kf.sample(T, initial_state=initial_state_mean)

# label half of the observations as missing
observations_missing = ma.array(
    observations_all, 
    mask=np.zeros(observations_all.shape)
)
for t in range(T):
    if t % 5 != 0:
        observations_missing[t] = ma.masked

# estimate state with filtering and smoothing
smoothed_states_all = kf.predict(observations_all)
smoothed_states_missing = kf.predict(observations_missing)

# draw estimates
pl.figure()
pl.hold(True)
lines_true = pl.plot(states, color='b')
lines_smooth_all = pl.plot(smoothed_states_all, color='r')
lines_smooth_missing = pl.plot(smoothed_states_missing, color='g')
pl.legend(
    (lines_true[0], lines_smooth_all[0], lines_smooth_missing[0]),
    ('true', 'all', 'missing'), loc='lower right'
)
pl.show()
