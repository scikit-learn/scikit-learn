'''
==================
Plot Kalman Simple
==================

This simple example shows how one may use the :ref:`sklearn.kalman`
module.
'''
import numpy as np
import matplotlib.pyplot as plt
from sklearn.kalman import KalmanFilter

# specify parameters
A = [[1, 0.1], [0, 1]]
b = [0, 0]
C = [1, -0.3]
d = 0
Q = np.eye(2)
R = 1
x_0 = [0, 0.2]
V_0 = [[1, 0.1], [-0.1, 1]]

# sample from model
rng = np.random.RandomState(0)
kf = KalmanFilter(A, C, Q, R, b, d, x_0, V_0, rng)
(x, z) = kf.sample(T=50)

# add a little noise and estimate model parameters
kf.Q += 0.1 * rng.multivariate_normal([0, 0], np.eye(2))
kf.R += 0.1 * rng.randn()
ll = kf.em(z, n_iter=5)[-1]
print 'log likelihood of observations for each EM iteration: {}'.format(ll)

# estimate state with filtering and smoothing
(x_filt, V_filt, loglik) = kf.filter(z)
(x_smooth, V_smooth, _) = kf.smooth(z)

# draw estimates
plt.figure()
plt.hold(True)
lines_true = plt.plot(x, color='b')
lines_filt = plt.plot(x_filt, color='r')
lines_smooth = plt.plot(x_smooth, color='g')
plt.legend((lines_true[0], lines_filt[0], lines_smooth[0]),
            ('true', 'filt', 'smooth'))
plt.show()
