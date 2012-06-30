import numpy as np
from numpy.testing import assert_almost_equal

from nose.tools import assert_true

from ..kalman import KalmanFilter
from sklearn.datasets import load_kalman_data

data = load_kalman_data()


def test_kalman_sampling():
    kf = KalmanFilter(
        data.transition_matrix,
        data.observation_matrix,
        data.transition_covariance,
        data.observation_covariance,
        data.transition_offsets,
        data.observation_offset,
        data.initial_state_mean,
        data.initial_state_covariance)

    (x, z) = kf.sample(100)
    assert_true(x.shape == (100, data.transition_matrix.shape[0]))
    assert_true(z.shape == (100, data.observation_matrix.shape[0]))


def test_kalman_filter_update():
    kf = KalmanFilter(
        data.transition_matrix,
        data.observation_matrix,
        data.transition_covariance,
        data.observation_covariance,
        data.transition_offsets,
        data.observation_offset,
        data.initial_state_mean,
        data.initial_state_covariance)

    # use Kalman Filter
    (x_filt, V_filt, ll) = kf.filter(X=data.data)

    # use online Kalman Filter
    T = data.data.shape[0]
    n_dim_state = data.transition_matrix.shape[0]
    x_filt2 = np.zeros((T, n_dim_state))
    V_filt2 = np.zeros((T, n_dim_state, n_dim_state))
    for t in range(T - 1):
        if t == 0:
            x_filt2[0] = data.initial_state_mean
            V_filt2[0] = data.initial_state_covariance
        (x_filt2[t + 1], V_filt2[t + 1], _) = kf.filter_update(
            x_filt2[t], V_filt2[t], data.data[t + 1], t=t
        )
    assert_true(np.all(x_filt == x_filt2))
    assert_true(np.all(V_filt == V_filt2))


def test_kalman_filter():
    kf = KalmanFilter(
        data.transition_matrix,
        data.observation_matrix,
        data.transition_covariance,
        data.observation_covariance,
        data.transition_offsets,
        data.observation_offset,
        data.initial_state_mean,
        data.initial_state_covariance)

    (x_filt, V_filt, ll) = kf.filter(X=data.data)
    for t in range(500):
        assert_true(np.linalg.norm(x_filt[t] - data.filtered_state_means[t]) < 1e-3)
        assert_true(np.linalg.norm(V_filt[t] - data.filtered_state_covariances[t]) < 1e-3)


def test_kalman_predict():
    kf = KalmanFilter(
        data.transition_matrix,
        data.observation_matrix,
        data.transition_covariance,
        data.observation_covariance,
        data.transition_offsets,
        data.observation_offset,
        data.initial_state_mean,
        data.initial_state_covariance)

    x_smooth = kf.predict(X=data.data)
    for t in reversed(range(501)):
        assert_true(np.linalg.norm(x_smooth[t] - data.smoothed_state_means[t]) < 1e-3)


def test_kalman_fit():
    # check against MATLAB dataset
    kf = KalmanFilter(
        data.transition_matrix,
        data.observation_matrix,
        data.initial_transition_covariance, 
        data.initial_observation_covariance,
        data.transition_offsets,
        data.observation_offset,
        data.initial_state_mean,
        data.initial_state_covariance,
        em_vars=['transition_covariance', 'observation_covariance'])

    scores = np.zeros(5)
    for i in range(len(scores)):
        scores[i] = np.sum(kf.filter(X=data.data)[-1])
        kf.fit(X=data.data, n_iter=1)

    assert_true(np.allclose(scores, data.loglikelihoods[:5]))

    # check that EM for all parameters is working
    kf.em_vars = 'all'
    T = 30
    for i in range(len(scores)):
        kf.fit(X=data.data[0:T], n_iter=1)
        scores[i] = np.sum(kf.filter(X=data.data[0:T])[-1])
    for i in range(len(scores) - 1):
        assert_true(scores[i] < scores[i + 1])
