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
kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q_0, R=data.R_0, b=data.b,
                  d=data.d, mu_0=data.x_0, sigma_0=data.V_0)

# Estimate mean and covariance of hidden state distribution iteratively.
T = data.data.shape[0]
n_dim_state = data.A.shape[0]
x_filt = np.zeros((T, n_dim_state))
V_filt = np.zeros((T, n_dim_state, n_dim_state))
for t in range(T - 1):
    if t == 0:
        x_filt[t] = data.x_0
        V_filt[t] = data.V_0
    (x_filt[t + 1], V_filt[t + 1], _) = kf.filter_update(
        x_filt[t], V_filt[t], data.data[t + 1], b=data.b[t], d=data.d
    )

# draw estimates
pl.figure()
pl.hold(True)
lines_true = pl.plot(data.target, color='b')
lines_filt = pl.plot(x_filt, color='r')
pl.legend((lines_true[0], lines_filt[0]),
            ('true',        'filt'))
pl.show()
