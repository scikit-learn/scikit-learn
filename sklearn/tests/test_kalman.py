import numpy as np

from ..kalman import KalmanFilter
from sklearn.datasets import load_kalman_data

data = load_kalman_data()


def test_kalman_sampling():
    kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q, R=data.R, b=data.b,
                      d=data.d, x_0=data.x_0, V_0=data.V_0)

    (x, z) = kf.sample(100)
    assert x.shape == (101, data.A.shape[0])
    assert z.shape == (100, data.C.shape[0])


def test_kalman_filter():
    kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q, R=data.R, b=data.b,
                      d=data.d, x_0=data.x_0, V_0=data.V_0)

    (x_filt, V_filt, ll) = kf.filter(Z=data.data)
    for t in range(500):
        assert np.linalg.norm(x_filt[t] - data.X_filt[t]) < 1e-5
        assert np.linalg.norm(V_filt[t] - data.V_filt[t]) < 1e-5


def test_kalman_smoother():
    kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q, R=data.R, b=data.b,
                      d=data.d, x_0=data.x_0, V_0=data.V_0)

    (x_smooth, V_smooth, ll) = kf.smooth(Z=data.data)
    for t in reversed(range(501)):
        assert np.linalg.norm(x_smooth[t] - data.X_smooth[t]) < 1e-5
        assert np.linalg.norm(V_smooth[t] - data.V_smooth[t]) < 1e-5


def test_kalman_em():
    kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q_0, R=data.R_0, b=data.b,
                      d=data.d, x_0=data.x_0, V_0=data.V_0)
    (Q, R, ll) = kf.em(Z=data.data, n_iter=5)
    assert np.allclose(ll, data.ll[:5])
