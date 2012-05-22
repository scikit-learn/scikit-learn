'''
===========================================
Using the Kalman Filter and Kalman Smoother
===========================================

This simple example shows how one may apply the Kalman Filter and Kalman
Smoother to some randomly generated data.
'''
import numpy as np
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

# sample from model
kf = KalmanFilter(A, C, Q, R, b, d, x_0, V_0, 
    random_state=np.random.RandomState(0))
(x, z) = kf.sample(T=50, x_0=x_0)

# estimate state with filtering and smoothing
x_filt = kf.filter(z)[0]
x_smooth = kf.predict(z)

# draw estimates
plt.figure()
plt.hold(True)
lines_true = plt.plot(x, color='b')
lines_filt = plt.plot(x_filt, color='r')
lines_smooth = plt.plot(x_smooth, color='g')
plt.legend((lines_true[0], lines_filt[0], lines_smooth[0]),
           (       'true',        'filt',        'smooth'),
           loc='lower right')
plt.show()
