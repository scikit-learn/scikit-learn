import pickle
from io import BytesIO
from StringIO import StringIO

import numpy as np
from numpy.testing import assert_array_almost_equal
from scipy import linalg
from nose.tools import assert_true

from sklearn.kalman import KalmanFilter
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
    (x_filt, V_filt) = kf.filter(X=data.data)

    # use online Kalman Filter
    n_timesteps = data.data.shape[0]
    n_dim_state = data.transition_matrix.shape[0]
    x_filt2 = np.zeros((n_timesteps, n_dim_state))
    V_filt2 = np.zeros((n_timesteps, n_dim_state, n_dim_state))
    for t in range(n_timesteps - 1):
        if t == 0:
            x_filt2[0] = data.initial_state_mean
            V_filt2[0] = data.initial_state_covariance
        (x_filt2[t + 1], V_filt2[t + 1]) = kf.filter_update(
            x_filt2[t], V_filt2[t], data.data[t + 1],
            transition_offset=data.transition_offsets[t]
        )
    assert_array_almost_equal(x_filt, x_filt2)
    assert_array_almost_equal(V_filt, V_filt2)


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

    (x_filt, V_filt) = kf.filter(X=data.data)
    assert_array_almost_equal(
        x_filt[:500],
        data.filtered_state_means[:500],
        decimal=3
    )
    assert_array_almost_equal(
        V_filt[:500],
        data.filtered_state_covariances[:500],
        decimal=3
    )


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
    assert_array_almost_equal(
        x_smooth[:501],
        data.smoothed_state_means[:501],
        decimal=3
    )


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
        scores[i] = kf.score(data.data)
        kf.fit(X=data.data, n_iter=1)

    assert_true(np.allclose(scores, data.loglikelihoods[:5]))

    # check that EM for all parameters is working
    kf.em_vars = 'all'
    n_timesteps = 30
    for i in range(len(scores)):
        kf.fit(X=data.data[0:n_timesteps], n_iter=1)
        scores[i] = kf.score(data.data[0:n_timesteps])
    for i in range(len(scores) - 1):
        assert_true(scores[i] < scores[i + 1])


def test_kalman_initialize_parameters():
    def check_dims(n_dim_state, n_dim_obs, kwargs):
        kf = KalmanFilter(**kwargs)
        (transition_matrices, transition_offsets, transition_covariance,
         observation_matrices, observation_offsets, observation_covariance,
         initial_state_mean, initial_state_covariance) = (
            kf._initialize_parameters()
        )
        assert_true(transition_matrices.shape == (n_dim_state, n_dim_state))
        assert_true(transition_offsets.shape == (n_dim_state,))
        assert_true(transition_covariance.shape == (n_dim_state, n_dim_state))
        assert_true(observation_matrices.shape == (n_dim_obs, n_dim_state))
        assert_true(observation_offsets.shape == (n_dim_obs,))
        assert_true(observation_covariance.shape == (n_dim_obs, n_dim_obs))
        assert_true(initial_state_mean.shape == (n_dim_state,))
        assert_true(
            initial_state_covariance.shape == (n_dim_state, n_dim_state)
        )

    check_dims(5, 1, {'transition_matrices': np.eye(5)})
    check_dims(1, 3, {'observation_offsets': np.zeros(3)})
    check_dims(2, 3, {'transition_covariance': np.eye(2),
                      'observation_offsets': np.zeros(3)})
    check_dims(3, 2, {'n_dim_state': 3, 'n_dim_obs': 2})
    check_dims(4, 1, {'initial_state_mean': np.zeros(4)})


def test_kalman_pickle():
    kf = KalmanFilter(
        data.transition_matrix,
        data.observation_matrix,
        data.transition_covariance,
        data.observation_covariance,
        data.transition_offsets,
        data.observation_offset,
        data.initial_state_mean,
        data.initial_state_covariance,
        em_vars='all')

    # train and get log likelihood
    X = data.data[0:10]
    kf = kf.fit(X, n_iter=5)
    score = kf.score(X)

    # pickle Kalman Filter
    store = StringIO()
    pickle.dump(kf, store)
    clf = pickle.load(StringIO(store.getvalue()))

    # check that parameters came out already
    np.testing.assert_almost_equal(score, kf.score(X))

    # store it as BytesIO as well
    store = BytesIO()
    pickle.dump(kf, store)
    kf = pickle.load(BytesIO(store.getvalue()))
