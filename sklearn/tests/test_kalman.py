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


def test_kalman_filter_update():
    kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q, R=data.R, b=data.b,
                      d=data.d, x_0=data.x_0, V_0=data.V_0)
    (x_filt, V_filt, ll) = kf.filter(Z=data.data)

    T = data.data.shape[0]
    n_dim_state = data.A.shape[0]
    x_filt2 = np.zeros((T + 1, n_dim_state))
    V_filt2 = np.zeros((T + 1, n_dim_state, n_dim_state))
    for t in range(T):
        if t == 0:
            x_filt2[0] = data.x_0
            V_filt2[0] = data.V_0
        (x_filt2[t+1], V_filt2[t+1], _) = kf.filter_update(x_filt[t],
            V_filt[t], data.data[t], t=t)
    assert np.all(x_filt == x_filt2)
    assert np.all(V_filt == V_filt2)


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
    # check against MATLAB dataset
    kf = KalmanFilter(A=data.A, C=data.C, Q=data.Q_0, R=data.R_0, b=data.b,
                      d=data.d, x_0=data.x_0, V_0=data.V_0, em_vars=['Q', 'R'])
    ll = kf.em(Z=data.data, n_iter=5)[-1]

    assert np.allclose(ll, data.ll[:5])

    # check that EM for all parameters is working
    kf.em_vars = 'all'
    ll = kf.em(Z=data.data[:40], n_iter=5)[-1]
    for i in range(len(ll)-1):
        assert ll[i] < ll[i+1]
