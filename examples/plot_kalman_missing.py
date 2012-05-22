'''
====================================================
Applying the Kalman Filter with Missing Observations
====================================================

This example shows how one may apply all of :mod:`sklearn.kalman`'s Kalman
Smooth, even with missing observations.
'''
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from sklearn.kalman import KalmanFilter

# specify parameters
random_state = np.random.RandomState(0)
A = [[1, 0.1], [0, 1]]
b = [-0.1, 0.1]
C = np.eye(2) + random_state.randn(2,2)*0.1
d = [1.0, -1.0]
Q = np.eye(2)
R = np.eye(2) + random_state.randn(2,2)*0.1
x_0 = [5, -5]
V_0 = [[1, 0.1], [-0.1, 1]]
T = 50

# sample from model
kf = KalmanFilter(A, C, Q, R, b, d, x_0, V_0, 
    random_state=np.random.RandomState(0))
(x, z_all) = kf.sample(T, x_0=x_0)

# label half of the observations as missing
z_missing = ma.array(z_all, mask=np.zeros(z_all.shape))
for t in range(T):
    if t % 5 != 0:
        z_missing.mask[t] = True

# estimate state with filtering and smoothing
x_smooth_all = kf.predict(z_all)
x_smooth_missing = kf.predict(z_missing)

# draw estimates
plt.figure()
plt.hold(True)
lines_true = plt.plot(x, color='b')
lines_smooth_all = plt.plot(x_smooth_all, color='r')
lines_smooth_missing = plt.plot(x_smooth_missing, color='g')
plt.legend(
    (lines_true[0], lines_smooth_all[0], lines_smooth_missing[0]),
    ('true', 'all', 'missing'), loc='lower right'
)
plt.show()

