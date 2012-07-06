'''
==================================
Kalman Filter tracking a sine wave
==================================

In this example, we generate a fake target trajectory using a sine wave.
Instead of observing those positions exactly, we observe the position plus some
random noise.  We then use a Kalman Filter to estimate the velocity of the
system as well.
'''
import numpy as np
import pylab as pl

from sklearn.kalman import KalmanFilter


rnd = np.random.RandomState(0)

# generate a noisy sine wave to act as our fake measurements
T = 100
x = np.linspace(0, 3 * np.pi, T)
observations = 20 * (np.sin(x) + 0.5 * rnd.randn(T))

# create a Kalman Filter by hinting at the size of the state and observation
# space.  If you already have good guesses for the initial parameters, put them
# in here.  The Kalman Filter will try to learn the values of all variables.
kf = KalmanFilter(
    transition_matrices=np.array([[1, 1], [0, 1]]),
    transition_covariance=np.eye(2) * 0.01,
    observation_covariance=10.0,
)

# You can use the Kalman Filter immediately without fitting, but its estimates
# may not be as good as if you fit first.
states_pred = kf.fit(observations).predict(observations)
print 'fitted model: %s' % (kf,)

# Plot lines for the measurements without noise, the estimated position of the
# target before fitting, and the estimated position after fitting.
pl.figure(figsize=(16, 6))
pl.hold(True)
obs_line = pl.plot(x, observations, linestyle='-', marker='x',
        color='b')
position_line = pl.plot(
    x, states_pred[:, 0], linestyle='-', marker='o', color='r'
)
velocity_line = pl.plot(
    x, states_pred[:, 1], linestyle='-', marker='o', color='g'
)
pl.legend(
    (obs_line[0], position_line[0], velocity_line[0]),
    ('true', 'position est.', 'velocity est.'),
    loc='lower right'
)
pl.xlim(xmax=x.max())
pl.show()
