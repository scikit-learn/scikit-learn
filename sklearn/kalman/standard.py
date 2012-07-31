"""
=====================================
Inference for Linear-Gaussian Systems
=====================================

This module implements the Kalman Filter, Kalman Smoother, and
EM Algorithm for Linear-Gaussian state space models
"""
import warnings

import numpy as np
from scipy import linalg

from ..base import BaseEstimator
from ..mixture import log_multivariate_normal_density
from ..utils import array1d, array2d, check_random_state

# Dimensionality of each Kalman Filter parameter for a single time step
DIM = {
    'transition_matrices': 2,
    'transition_offsets': 1,
    'observation_matrices': 2,
    'observation_offsets': 1,
    'transition_covariance': 2,
    'observation_covariance': 2,
    'initial_state_mean': 1,
    'initial_state_covariance': 2,
}


def _arg_or_default(arg, default, dim, name):
    if arg is None:
        result = default
    else:
        result = arg
    if len(result.shape) > dim:
        raise ValueError(
            ('%s is not constant for all time.'
             + '  You must specify it manually.') % (name,)
        )
    return result


def _determine_dimensionality(variables, default):
    """Derive the dimensionality of the state space

    Parameters
    ----------
    variables : list of ({None, array}, conversion function, index)
        variables, functions to convert them to arrays, and indices in those
        arrays to derive dimensionality from.
    default : {None, int}
        default dimensionality to return if variables is empty

    Returns
    -------
    dim : int
        dimensionality of state space as derived from variables or default.
    """
    # gather possible values based on the variables
    candidates = []
    for (v, converter, idx) in variables:
        if v is not None:
            v = converter(v)
            candidates.append(v.shape[idx])

    # also use the manually specified default
    if default is not None:
        candidates.append(default)

    # ensure consistency of all derived values
    if len(candidates) == 0:
        return 1
    else:
        if not np.all(np.array(candidates) == candidates[0]):
            raise ValueError(
                "The shape of all " +
                "parameters is not consistent.  " +
                "Please re-check their values."
            )
        return candidates[0]


def _last_dims(X, t, ndims=2):
    """Extract the final dimensions of `X`

    Extract the final `ndim` dimensions at index `t` if `X` has >= `ndim` + 1
    dimensions, otherwise return `X`.

    Parameters
    ----------
    X : array with at least dimension `ndims`
    t : int
        index to use for the `ndims` + 1th dimension
    ndims : int, optional
        number of dimensions in the array desired

    Returns
    -------
    Y : array with dimension `ndims`
        the final `ndims` dimensions indexed by `t`
    """
    X = np.asarray(X)
    if len(X.shape) == ndims + 1:
        return X[t]
    elif len(X.shape) == ndims:
        return X
    else:
        raise ValueError(("X only has %d dimensions when %d" +
                " or more are required") % (len(X.shape), ndims))


def _loglikelihoods(observation_matrices, observation_offsets,
                    observation_covariance, predicted_state_means,
                    predicted_state_covariances, observations):
    """Calculate log likelihood of all observations

    Parameters
    ----------
    observation_matrices : [n_timesteps, n_dim_obs, n_dim_obs] or [n_dim_obs,
    n_dim_state] array
        observation matrices for t in [0...n_timesteps-1]
    observation_offsets : [n_timesteps, n_dim_obs] or [n_dim_obs] array
        offsets for observations for t = [0...n_timesteps-1]
    observation_covariance : [n_dim_obs, n_dim_obs] array
        covariance matrix for all observations
    predicted_state_means : [n_timesteps, n_dim_state] array
        mean of state at time t given observations from times
        [0...t-1] for t in [0...n_timesteps-1]
    predicted_state_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [0...t-1] for t in [0...n_timesteps-1]
    observations : [n_dim_obs] array
        All observations.  If `observations[t]` is a masked array and any of
        its values are masked, the observation will be ignored.

    Returns
    -------
    loglikelihoods: [n_timesteps] array
        `loglikelihoods[t]` is the probability density of the observation
        generated at time step t
    """
    n_timesteps = observations.shape[0]
    loglikelihoods = np.zeros(n_timesteps)
    for t in range(n_timesteps):
        observation = observations[t]
        if not np.any(np.ma.getmask(observation)):
            observation_matrix = _last_dims(observation_matrices, t)
            observation_offset = _last_dims(observation_offsets, t, ndims=1)
            predicted_state_mean = _last_dims(
                predicted_state_means, t, ndims=1
            )
            predicted_state_covariance = _last_dims(
                predicted_state_covariances, t
            )

            predicted_observation_mean = (
                np.dot(observation_matrix,
                       predicted_state_mean)
                + observation_offset
            )
            predicted_observation_covariance = (
                np.dot(observation_matrix,
                       np.dot(predicted_state_covariance,
                              observation_matrix.T))
                + observation_covariance
            )
            loglikelihoods[t] = log_multivariate_normal_density(
                observation[np.newaxis, :],
                predicted_observation_mean[np.newaxis, :],
                predicted_observation_covariance[np.newaxis, :, :],
                'full'
            )
    return loglikelihoods


def _filter_predict(transition_matrix, transition_covariance,
                    transition_offset, current_state_mean,
                    current_state_covariance):
    r"""Calculate the mean and covariance of :math:`P(x_{t+1} | z_{0:t})`

    Using the mean and covariance of :math:`P(x_t | z_{0:t})`, calculate the
    mean and covariance of :math:`P(x_{t+1} | z_{0:t})`.

    Parameters
    ----------
    transition_matrix : [n_dim_state, n_dim_state} array
        state transition matrix from time t to t+1
    transition_covariance : [n_dim_state, n_dim_state] array
        covariance matrix for state transition from time t to t+1
    transition_offset : [n_dim_state] array
        offset for state transition from time t to t+1
    current_state_mean: [n_dim_state] array
        mean of state at time t given observations from times
        [0...t]
    current_state_covariance: [n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [0...t]

    Returns
    -------
    predicted_state_mean : [n_dim_state] array
        mean of state at time t+1 given observations from times [0...t]
    predicted_state_covariance : [n_dim_state, n_dim_state] array
        covariance of state at time t+1 given observations from times
        [0...t]
    """
    predicted_state_mean = (
        np.dot(transition_matrix, current_state_mean)
        + transition_offset
    )
    predicted_state_covariance = (
        np.dot(transition_matrix,
               np.dot(current_state_covariance,
                      transition_matrix.T))
        + transition_covariance
    )

    return (predicted_state_mean, predicted_state_covariance)


def _filter_correct(observation_matrix, observation_covariance,
                    observation_offset, predicted_state_mean,
                    predicted_state_covariance, observation):
    r"""Correct a predicted state with a Kalman Filter update

    Incorporate observation `observation` from time `t` to turn
    :math:`P(x_t | z_{0:t-1})` into :math:`P(x_t | z_{0:t})`

    Parameters
    ----------
    observation_matrix : [n_dim_obs, n_dim_state] array
        observation matrix for time t
    observation_covariance : [n_dim_obs, n_dim_obs] array
        covariance matrix for observation at time t
    observation_offset : [n_dim_obs] array
        offset for observation at time t
    predicted_state_mean : [n_dim_state] array
        mean of state at time t given observations from times
        [0...t-1]
    predicted_state_covariance : [n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [0...t-1]
    observation : [n_dim_obs] array
        observation at time t.  If `observation` is a masked array and any of
        its values are masked, the observation will be ignored.

    Returns
    -------
    kalman_gain : [n_dim_state, n_dim_obs] array
        Kalman gain matrix for time t
    corrected_state_mean : [n_dim_state] array
        mean of state at time t given observations from times
        [0...t]
    corrected_state_covariance : [n_dim_state, n_dim_state] array
        covariance of state at time t given observations from times
        [0...t]
    """
    if not np.any(np.ma.getmask(observation)):
        predicted_observation_mean = (
            np.dot(observation_matrix,
                   predicted_state_mean)
            + observation_offset
        )
        predicted_observation_covariance = (
            np.dot(observation_matrix,
                   np.dot(predicted_state_covariance,
                          observation_matrix.T))
            + observation_covariance
        )

        kalman_gain = (
            np.dot(predicted_state_covariance,
                   np.dot(observation_matrix.T,
                          linalg.pinv(predicted_observation_covariance)))
        )

        corrected_state_mean = (
            predicted_state_mean
            + np.dot(kalman_gain, observation - predicted_observation_mean)
        )
        corrected_state_covariance = (
            predicted_state_covariance
            - np.dot(kalman_gain,
                     np.dot(observation_matrix,
                            predicted_state_covariance))
        )
    else:
        n_dim_state = predicted_state_covariance.shape[0]
        n_dim_obs = observation_matrix.shape[0]
        kalman_gain = np.zeros((n_dim_state, n_dim_obs))

        corrected_state_mean = predicted_state_mean
        corrected_state_covariance = predicted_state_covariance

    return (kalman_gain, corrected_state_mean,
            corrected_state_covariance)


def _filter(transition_matrices, observation_matrices, transition_covariance,
            observation_covariance, transition_offsets, observation_offsets,
            initial_state_mean, initial_state_covariance, observations):
    """Apply the Kalman Filter

    Calculate posterior distribution over hidden states given observations up
    to and including the current time step.

    Parameters
    ----------
    transition_matrices : [n_timesteps-1,n_dim_state,n_dim_state] or
    [n_dim_state,n_dim_state] array-like
        state transition matrices
    observation_matrices : [n_timesteps, n_dim_obs, n_dim_obs] or [n_dim_obs, \
    n_dim_obs] array-like
        observation matrix
    transition_covariance : [n_timesteps-1,n_dim_state,n_dim_state] or
    [n_dim_state,n_dim_state] array-like
        state transition covariance matrix
    observation_covariance : [n_timesteps, n_dim_obs, n_dim_obs] or [n_dim_obs,
    n_dim_obs] array-like
        observation covariance matrix
    transition_offsets : [n_timesteps-1, n_dim_state] or [n_dim_state] \
    array-like
        state offset
    observation_offsets : [n_timesteps, n_dim_obs] or [n_dim_obs] array-like
        observations for times [0...n_timesteps-1]
    initial_state_mean : [n_dim_state] array-like
        mean of initial state distribution
    initial_state_covariance : [n_dim_state, n_dim_state] array-like
        covariance of initial state distribution
    observations : [n_timesteps, n_dim_obs] array
        observations from times [0...n_timesteps-1].  If `observations` is a
        masked array and any of `observations[t]` is masked, then
        `observations[t]` will be treated as a missing observation.

    Returns
    -------
    predicted_state_means : [n_timesteps, n_dim_state] array
        `predicted_state_means[t]` = mean of hidden state at time t given
        observations from times [0...t-1]
    predicted_state_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        `predicted_state_covariances[t]` = covariance of hidden state at time t
        given observations from times [0...t-1]
    kalman_gains : [n_timesteps, n_dim_state] array
        `kalman_gains[t]` = Kalman gain matrix for time t
    filtered_state_means : [n_timesteps, n_dim_state] array
        `filtered_state_means[t]` = mean of hidden state at time t given
        observations from times [0...t]
    filtered_state_covariances : [n_timesteps, n_dim_state] array
        `filtered_state_covariances[t]` = covariance of hidden state at time t
        given observations from times [0...t]
    """
    n_timesteps = observations.shape[0]
    n_dim_state = len(initial_state_mean)
    n_dim_obs = observations.shape[1]

    predicted_state_means = np.zeros((n_timesteps, n_dim_state))
    predicted_state_covariances = np.zeros(
        (n_timesteps, n_dim_state, n_dim_state)
    )
    kalman_gains = np.zeros((n_timesteps, n_dim_state, n_dim_obs))
    filtered_state_means = np.zeros((n_timesteps, n_dim_state))
    filtered_state_covariances = np.zeros(
        (n_timesteps, n_dim_state, n_dim_state)
    )

    for t in range(n_timesteps):
        if t == 0:
            predicted_state_means[t] = initial_state_mean
            predicted_state_covariances[t] = initial_state_covariance
        else:
            transition_matrix = _last_dims(transition_matrices, t - 1)
            transition_covariance = _last_dims(transition_covariance, t - 1)
            transition_offset = _last_dims(transition_offsets, t - 1, ndims=1)
            predicted_state_means[t], predicted_state_covariances[t] = (
                _filter_predict(
                    transition_matrix,
                    transition_covariance,
                    transition_offset,
                    filtered_state_means[t - 1],
                    filtered_state_covariances[t - 1]
                )
            )

        observation_matrix = _last_dims(observation_matrices, t)
        observation_covariance = _last_dims(observation_covariance, t)
        observation_offset = _last_dims(observation_offsets, t, ndims=1)
        (kalman_gains[t], filtered_state_means[t],
         filtered_state_covariances[t]) = (
            _filter_correct(observation_matrix,
                observation_covariance,
                observation_offset,
                predicted_state_means[t],
                predicted_state_covariances[t],
                observations[t]
            )
        )

    return (predicted_state_means, predicted_state_covariances,
            kalman_gains, filtered_state_means,
            filtered_state_covariances)


def _smooth_update(transition_matrix, filtered_state_mean,
                   filtered_state_covariance, predicted_state_mean,
                   predicted_state_covariance, next_smoothed_state_mean,
                   next_smoothed_state_covariance):
    r"""Correct a predicted state with a Kalman Smoother update

    Calculates posterior distribution of the hidden state at time `t` given the
    observations all observations via Kalman Smoothing.

    Parameters
    ----------
    transition_matrix : [n_dim_state, n_dim_state] array
        state transition matrix from time t to t+1
    filtered_state_mean : [n_dim_state] array
        mean of filtered state at time t given observations from
        times [0...t]
    filtered_state_covariance : [n_dim_state, n_dim_state] array
        covariance of filtered state at time t given observations from
        times [0...t]
    predicted_state_mean : [n_dim_state] array
        mean of filtered state at time t+1 given observations from
        times [0...t]
    predicted_state_covariance : [n_dim_state, n_dim_state] array
        covariance of filtered state at time t+1 given observations from
        times [0...t]
    next_smoothed_state_mean : [n_dim_state] array
        mean of smoothed state at time t+1 given observations from
        times [0...n_timesteps-1]
    next_smoothed_state_covariance : [n_dim_state, n_dim_state] array
        covariance of smoothed state at time t+1 given observations from
        times [0...n_timesteps-1]

    Returns
    -------
    smoothed_state_mean : [n_dim_state] array
        mean of smoothed state at time t given observations from times
        [0...n_timesteps-1]
    smoothed_state_covariance : [n_dim_state, n_dim_state] array
        covariance of smoothed state at time t given observations from
        times [0...n_timesteps-1]
    kalman_smoothing_gain : [n_dim_state, n_dim_state] array
        correction matrix for Kalman Smoothing at time t
    """
    kalman_smoothing_gain = (
        np.dot(filtered_state_covariance,
               np.dot(transition_matrix.T,
                      linalg.pinv(predicted_state_covariance)))
    )

    smoothed_state_mean = (
        filtered_state_mean
        + np.dot(kalman_smoothing_gain,
                 next_smoothed_state_mean - predicted_state_mean)
    )
    smoothed_state_covariance = (
        filtered_state_covariance
        + np.dot(kalman_smoothing_gain,
                 np.dot(
                    (next_smoothed_state_covariance
                        - predicted_state_covariance),
                    kalman_smoothing_gain.T
                 ))
    )

    return (smoothed_state_mean, smoothed_state_covariance,
            kalman_smoothing_gain)


def _smooth(transition_matrices, filtered_state_means,
            filtered_state_covariances, predicted_state_means,
            predicted_state_covariances):
    """Apply the Kalman Smoother

    Estimate the hidden state at time for each time step given all
    observations.

    Parameters
    ----------
    transition_matrices : [n_timesteps-1, n_dim_state, n_dim_state] or \
    [n_dim_state, n_dim_state] array
        `transition_matrices[t]` = transition matrix from time t to t+1
    filtered_state_means : [n_timesteps, n_dim_state] array
        `filtered_state_means[t]` = mean state estimate for time t given
        observations from times [0...t]
    filtered_state_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        `filtered_state_covariances[t]` = covariance of state estimate for time
        t given observations from times [0...t]
    predicted_state_means : [n_timesteps, n_dim_state] array
        `predicted_state_means[t]` = mean state estimate for time t given
        observations from times [0...t-1]
    predicted_state_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        `predicted_state_covariances[t]` = covariance of state estimate for
        time t given observations from times [0...t-1]

    Returns
    -------
    smoothed_state_means : [n_timesteps, n_dim_state]
        mean of hidden state distributions for times [0...n_timesteps-1] given
        all observations
    smoothed_state_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        covariance matrix of hidden state distributions for times
        [0...n_timesteps-1] given all observations
    kalman_smoothing_gains : [n_timesteps-1, n_dim_state, n_dim_state] array
        Kalman Smoothing correction matrices for times [0...n_timesteps-2]
    """
    n_timesteps, n_dim_state = filtered_state_means.shape

    smoothed_state_means = np.zeros((n_timesteps, n_dim_state))
    smoothed_state_covariances = np.zeros((n_timesteps, n_dim_state,
                                           n_dim_state))
    kalman_smoothing_gains = np.zeros((n_timesteps - 1, n_dim_state,
                                       n_dim_state))

    smoothed_state_means[-1] = filtered_state_means[-1]
    smoothed_state_covariances[-1] = filtered_state_covariances[-1]

    for t in reversed(range(n_timesteps - 1)):
        transition_matrix = _last_dims(transition_matrices, t)
        (smoothed_state_means[t], smoothed_state_covariances[t],
         kalman_smoothing_gains[t]) = (
            _smooth_update(
                transition_matrix,
                filtered_state_means[t],
                filtered_state_covariances[t],
                predicted_state_means[t + 1],
                predicted_state_covariances[t + 1],
                smoothed_state_means[t + 1],
                smoothed_state_covariances[t + 1]
            )
        )
    return (smoothed_state_means, smoothed_state_covariances,
            kalman_smoothing_gains)


def _smooth_pair(smoothed_state_covariances, kalman_smoothing_gain):
    r"""Calculate pairwise covariance between hidden states

    Calculate covariance between hidden states at :math:`t` and :math:`t-1` for
    all time step pairs

    Parameters
    ----------
    smoothed_state_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        covariance of hidden state given all observations
    kalman_smoothing_gain : [n_timesteps-1, n_dim_state, n_dim_state]
        Correction matrices from Kalman Smoothing

    Returns
    -------
    pairwise_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        Covariance between hidden states at times t and t-1 for t =
        [1...n_timesteps-1].  Time 0 is ignored.
    """
    n_timesteps, n_dim_state, _ = smoothed_state_covariances.shape
    pairwise_covariances = np.zeros((n_timesteps, n_dim_state, n_dim_state))
    for t in range(1, n_timesteps):
        pairwise_covariances[t] = (
            np.dot(smoothed_state_covariances[t],
                   kalman_smoothing_gain[t - 1].T)
        )
    return pairwise_covariances


def _em(observations, transition_offsets, observation_offsets,
        smoothed_state_means, smoothed_state_covariances, pairwise_covariances,
        given={}):
    """Apply the EM Algorithm to the Linear-Gaussian model

    Estimate Linear-Gaussian model parameters by maximizing the expected log
    likelihood of all observations.

    Parameters
    ----------
    observations : [n_timesteps, n_dim_obs] array
        observations for times [0...n_timesteps-1].  If observations is a
        masked array and any of observations[t] is masked, then it will be
        treated as a missing observation.
    transition_offsets : [n_dim_state] or [n_timesteps-1, n_dim_state] array
        transition offset
    observation_offsets : [n_dim_obs] or [n_timesteps, n_dim_obs] array
        observation offsets
    smoothed_state_means : [n_timesteps, n_dim_state] array
        smoothed_state_means[t] = mean of state at time t given all
        observations
    smoothed_state_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        smoothed_state_covariances[t] = covariance of state at time t given all
        observations
    pairwise_covariances : [n_timesteps, n_dim_state, n_dim_state] array
        pairwise_covariances[t] = covariance between states at times t and
        t-1 given all observations.  pairwise_covariances[0] is ignored.
    given: dict
        if one of the variables EM is capable of predicting is in given, then
        that value will be used and EM will not attempt to estimate it.  e.g.,
        if 'observation_matrix' is in given, observation_matrix will not be
        estimated and given['observation_matrix'] will be returned in its
        place.

    Returns
    -------
    transition_matrix : [n_dim_state, n_dim_state] array
        estimated transition matrix
    observation_matrix : [n_dim_obs, n_dim_state] array
        estimated observation matrix
    transition_offsets : [n_dim_state] array
        estimated transition offset
    observation_offsets : [n_dim_obs] array
        estimated observation offset
    transition_covariance : [n_dim_state, n_dim_state] array
        estimated covariance matrix for state transitions
    observation_covariance : [n_dim_obs, n_dim_obs] array
        estimated covariance matrix for observations
    initial_state_mean : [n_dim_state] array
        estimated mean of initial state distribution
    initial_state_covariance : [n_dim_state] array
        estimated covariance of initial state distribution
    """
    if 'observation_matrices' in given:
        observation_matrix = given['observation_matrices']
    else:
        observation_matrix = _em_observation_matrix(
            observations, observation_offsets,
            smoothed_state_means, smoothed_state_covariances
        )

    if 'observation_covariance' in given:
        observation_covariance = given['observation_covariance']
    else:
        observation_covariance = _em_observation_covariance(
            observations, observation_offsets,
            observation_matrix, smoothed_state_means,
            smoothed_state_covariances
        )

    if 'transition_matrices' in given:
        transition_matrix = given['transition_matrices']
    else:
        transition_matrix = _em_transition_matrix(
            transition_offsets, smoothed_state_means,
            smoothed_state_covariances, pairwise_covariances
        )

    if 'transition_covariance' in given:
        transition_covariance = given['transition_covariance']
    else:
        transition_covariance = _em_transition_covariance(
            transition_matrix, transition_offsets,
            smoothed_state_means, smoothed_state_covariances,
            pairwise_covariances
        )

    if 'initial_state_mean' in given:
        initial_state_mean = given['initial_state_mean']
    else:
        initial_state_mean = _em_initial_state_mean(smoothed_state_means)

    if 'initial_state_covariance' in given:
        initial_state_covariance = given['initial_state_covariance']
    else:
        initial_state_covariance = _em_initial_state_covariance(
            initial_state_mean, smoothed_state_means,
            smoothed_state_covariances
        )

    if 'transition_offsets' in given:
        transition_offset = given['transition_offsets']
    else:
        transition_offset = _em_transition_offset(
            transition_matrix,
            smoothed_state_means
        )

    if 'observation_offsets' in given:
        observation_offset = given['observation_offsets']
    else:
        observation_offset = _em_observation_offset(
            observation_matrix, smoothed_state_means,
            observations
        )

    return (transition_matrix, observation_matrix, transition_offset,
            observation_offset, transition_covariance,
            observation_covariance, initial_state_mean,
            initial_state_covariance)


def _em_observation_matrix(observations, observation_offsets,
                          smoothed_state_means, smoothed_state_covariances):
    r"""Apply the EM algorithm to parameter `observation_matrix`

    Maximize expected log likelihood of observations with respect to the
    observation matrix `observation_matrix`.

    .. math::

        C &= ( \sum_{t=0}^{T-1} (z_t - d_t) \mathbb{E}[x_t]^T )
             ( \sum_{t=0}^{T-1} \mathbb{E}[x_t x_t^T] )^-1

    """
    _, n_dim_state = smoothed_state_means.shape
    n_timesteps, n_dim_obs = observations.shape
    res1 = np.zeros((n_dim_obs, n_dim_state))
    res2 = np.zeros((n_dim_state, n_dim_state))
    for t in range(n_timesteps):
        if not np.any(np.ma.getmask(observations[t])):
            observation_offset = _last_dims(observation_offsets, t, ndims=1)
            res1 += np.outer(observations[t] - observation_offset,
                             smoothed_state_means[t])
            res2 += (
                smoothed_state_covariances[t]
                + np.outer(smoothed_state_means[t], smoothed_state_means[t])
            )
    return np.dot(res1, linalg.pinv(res2))


def _em_observation_covariance(observations, observation_offsets,
                              transition_matrices, smoothed_state_means,
                              smoothed_state_covariances):
    r"""Apply the EM algorithm to parameter `observation_covariance`

    Maximize expected log likelihood of observations with respect to the
    observation covariance matrix `observation_covariance`.

    .. math::

        R &= \frac{1}{T} \sum_{t=0}^{T-1}
                [z_t - C_t \mathbb{E}[x_t] - b_t]
                    [z_t - C_t \mathbb{E}[x_t] - b_t]^T
                + C_t Var(x_t) C_t^T
    """
    _, n_dim_state = smoothed_state_means.shape
    n_timesteps, n_dim_obs = observations.shape
    res = np.zeros((n_dim_obs, n_dim_obs))
    n_obs = 0
    for t in range(n_timesteps):
        if not np.any(np.ma.getmask(observations[t])):
            transition_matrix = _last_dims(transition_matrices, t)
            transition_offset = _last_dims(observation_offsets, t, ndims=1)
            err = (
                observations[t]
                - np.dot(transition_matrix, smoothed_state_means[t])
                - transition_offset
            )
            res += (
                np.outer(err, err)
                + np.dot(transition_matrix,
                         np.dot(smoothed_state_covariances[t],
                                transition_matrix.T))
            )
            n_obs += 1
    if n_obs > 0:
        return (1.0 / n_obs) * res
    else:
        return res


def _em_transition_matrix(transition_offsets, smoothed_state_means,
                          smoothed_state_covariances, pairwise_covariances):
    r"""Apply the EM algorithm to parameter `transition_matrix`

    Maximize expected log likelihood of observations with respect to the state
    transition matrix `transition_matrix`.

    .. math::

        A &= ( \sum_{t=1}^{T-1} \mathbb{E}[x_t x_{t-1}^{T}]
                - b_{t-1} \mathbb{E}[x_{t-1}]^T )
             ( \sum_{t=1}^{T-1} \mathbb{E}[x_{t-1} x_{t-1}^T] )^{-1}
    """
    n_timesteps, n_dim_state, _ = smoothed_state_covariances.shape
    res1 = np.zeros((n_dim_state, n_dim_state))
    res2 = np.zeros((n_dim_state, n_dim_state))
    for t in range(1, n_timesteps):
        transition_offset = _last_dims(transition_offsets, t - 1, ndims=1)
        res1 += (
            pairwise_covariances[t]
            + np.outer(smoothed_state_means[t],
                       smoothed_state_means[t - 1])
            - np.outer(transition_offset, smoothed_state_means[t - 1])
        )
        res2 += (
            smoothed_state_covariances[t - 1]
            + np.outer(smoothed_state_means[t - 1],
                       smoothed_state_means[t - 1])
        )
    return np.dot(res1, linalg.pinv(res2))


def _em_transition_covariance(transition_matrices, transition_offsets,
                              smoothed_state_means, smoothed_state_covariances,
                              pairwise_covariances):
    r"""Apply the EM algorithm to parameter `transition_covariance`

    Maximize expected log likelihood of observations with respect to the
    transition covariance matrix `transition_covariance`.

    .. math::

        Q &= \frac{1}{T-1} \sum_{t=0}^{T-2}
                (\mathbb{E}[x_{t+1}] - A_t \mathbb{E}[x_t] - b_t)
                    (\mathbb{E}[x_{t+1}] - A_t \mathbb{E}[x_t] - b_t)^T
                + A_t Var(x_t) A_t^T + Var(x_{t+1})
                - Cov(x_{t+1}, x_t) A_t^T - A_t Cov(x_t, x_{t+1})
    """
    n_timesteps, n_dim_state, _ = smoothed_state_covariances.shape
    res = np.zeros((n_dim_state, n_dim_state))
    for t in range(n_timesteps - 1):
        transition_matrix = _last_dims(transition_matrices, t)
        transition_offset = _last_dims(transition_offsets, t, ndims=1)
        err = (
            smoothed_state_means[t + 1]
            - np.dot(transition_matrix, smoothed_state_means[t])
            - transition_offset
        )
        Vt1t_A = (
            np.dot(pairwise_covariances[t + 1],
                   transition_matrix.T)
        )
        res += (
            np.outer(err, err)
            + np.dot(transition_matrix,
                     np.dot(smoothed_state_covariances[t],
                            transition_matrix.T))
            + smoothed_state_covariances[t + 1]
            - Vt1t_A - Vt1t_A.T
        )

    return (1.0 / (n_timesteps - 1)) * res


def _em_initial_state_mean(smoothed_state_means):
    r"""Apply the EM algorithm to parameter `initial_state_mean`

    Maximize expected log likelihood of observations with respect to the
    initial state distribution mean `initial_state_mean`.

    .. math::

        \mu_0 = \mathbb{E}[x_0]
    """

    return smoothed_state_means[0]


def _em_initial_state_covariance(initial_state_mean, smoothed_state_means,
                                 smoothed_state_covariances):
    r"""Apply the EM algorithm to parameter `initial_state_covariance`

    Maximize expected log likelihood of observations with respect to the
    covariance of the initial state distribution `initial_state_covariance`.

    .. math::

        \sigma_0 = \mathbb{E}[x_0, x_0^T] - mu_0 \mathbb{E}[x_0]^T
                   - \mathbb{E}[x_0] mu_0^T + mu_0 mu_0^T
    """
    x0 = smoothed_state_means[0]
    x0_x0 = smoothed_state_covariances[0] + np.outer(x0, x0)
    return (
        x0_x0
        - np.outer(initial_state_mean, x0)
        - np.outer(x0, initial_state_mean)
        + np.outer(initial_state_mean, initial_state_mean)
    )


def _em_transition_offset(transition_matrices, smoothed_state_means):
    r"""Apply the EM algorithm to parameter `transition_offset`

    Maximize expected log likelihood of observations with respect to the
    state transition offset `transition_offset`.

    .. math::

        b = \frac{1}{T-1} \sum_{t=1}^{T-1}
                \mathbb{E}[x_t] - A_{t-1} \mathbb{E}[x_{t-1}]
    """
    n_timesteps, n_dim_state = smoothed_state_means.shape
    transition_offset = np.zeros(n_dim_state)
    for t in range(1, n_timesteps):
        transition_matrix = _last_dims(transition_matrices, t - 1)
        transition_offset += (
            smoothed_state_means[t]
            - np.dot(transition_matrix, smoothed_state_means[t - 1])
        )
    if n_timesteps > 1:
        return (1.0 / (n_timesteps - 1)) * transition_offset
    else:
        return np.zeros(n_dim_state)


def _em_observation_offset(observation_matrices, smoothed_state_means,
                           observations):
    r"""Apply the EM algorithm to parameter `observation_offset`

    Maximize expected log likelihood of observations with respect to the
    observation offset `observation_offset`.

    .. math::

        d = \frac{1}{T} \sum_{t=0}^{T-1} z_t - C_{t} \mathbb{E}[x_{t}]
    """
    n_timesteps, n_dim_obs = observations.shape
    observation_offset = np.zeros(n_dim_obs)
    n_obs = 0
    for t in range(n_timesteps):
        if not np.any(np.ma.getmask(observations[t])):
            observation_matrix = _last_dims(observation_matrices, t)
            observation_offset += (
                observations[t]
                - np.dot(observation_matrix, smoothed_state_means[t])
            )
            n_obs += 1
    if n_obs > 0:
        return (1.0 / n_obs) * observation_offset
    else:
        return observation_offset


class KalmanFilter(BaseEstimator):
    """Implements the Kalman Filter, Kalman Smoother, and EM algorithm.

    This class implements the Kalman Filter, Kalman Smoother, and EM Algorithm
    for a Linear Gaussian model specified by,

    .. math::

        x_{t+1}   &= A_{t} x_{t} + b_{t} + \\text{Normal}(0, Q_{t}) \\\\
        z_{t}     &= C_{t} x_{t} + d_{t} + \\text{Normal}(0, R_{t})

    The Kalman Filter is an algorithm designed to estimate
    :math:`P(x_t | z_{0:t})`.  As all state transitions and observations are
    linear with Gaussian distributed noise, these distributions can be
    represented exactly as Gaussian distributions with mean `mu_filt[t]` and
    covariances `sigma_filt[t]`.

    Similarly, the Kalman Smoother is an algorithm designed to estimate
    :math:`P(x_t | z_{0:T-1})`.

    The EM algorithm aims to find for
    :math:`\\theta = (A, b, C, d, Q, R, \\mu_0, \\sigma_0)`

    .. math::

        \\max_{\\theta} P(z_{0:T-1}; \\theta)

    If we define :math:`L(x_{0:t},\\theta) = \\log P(z_{0:T-1}, x_{0:T-1};
    \\theta)`, then the EM algorithm works by iteratively finding,

    .. math::

        P(x_{0:T-1} | z_{0:T-1}, \\theta_i)

    then by maximizing,

    .. math::

        \\theta_{i+1} = \\arg\\max_{\\theta}
            \\mathbb{E}_{x_{0:T-1}} [
                L(x_{0:t}, \\theta)| z_{0:T-1}, \\theta_i
            ]

    Parameters
    ----------
    transition_matrices : [n_timesteps-1, n_dim_state, n_dim_state] or \
    [n_dim_state,n_dim_state] array-like
        Also known as :math:`A`.  state transition matrix between times t and
        t+1 for t in [0...n_timesteps-2]
    observation_matrices : [n_timesteps, n_dim_obs, n_dim_obs] or [n_dim_obs, \
    n_dim_obs] array-like
        Also known as :math:`C`.  observation matrix for times
        [0...n_timesteps-1]
    transition_covariance : [n_dim_state, n_dim_state] array-like
        Also known as :math:`Q`.  state transition covariance matrix for times
        [0...n_timesteps-2]
    observation_covariance : [n_dim_obs, n_dim_obs] array-like
        Also known as :math:`R`.  observation covariance matrix for times
        [0...n_timesteps-1]
    transition_offsets : [n_timesteps-1, n_dim_state] or [n_dim_state] \
    array-like
        Also known as :math:`b`.  state offsets for times [0...n_timesteps-2]
    observation_offsets : [n_timesteps, n_dim_obs] or [n_dim_obs] array-like
        Also known as :math:`d`.  observation offset for times
        [0...n_timesteps-1]
    initial_state_mean : [n_dim_state] array-like
        Also known as :math:`\mu_0`mean of initial state distribution
    initial_state_covariance : [n_dim_state, n_dim_state] array-like
        Also known as :math:`\sigma_0`.  covariance of initial state
        distribution
    random_state : optional, numpy random state
        random number generator used in sampling
    em_vars : optional, subset of ['transition_matrices', \
    'observation_matrices', 'transition_offsets', 'observation_offsets', \
    'transition_covariance', 'observation_covariance', 'initial_state_mean', \
    'initial_state_covariance'] or 'all'
        if `em_vars` is an iterable of strings only variables in `em_vars`
        will be estimated using EM.  if `em_vars` == 'all', then all
        variables will be estimated.
    n_dim_state: optional, integer
        the dimensionality of the state space. Only meaningful when you do not
        specify initial values for `transition_matrices`, `transition_offsets`,
        `transition_covariance`, `initial_state_mean`, or
        `initial_state_covariance`.
    n_dim_obs: optional, integer
        the dimensionality of the observation space. Only meaningful when you
        do not specify initial values for `observation_matrices`,
        `observation_offsets`, or `observation_covariance`.
    """
    def __init__(self, transition_matrices=None, observation_matrices=None,
            transition_covariance=None, observation_covariance=None,
            transition_offsets=None, observation_offsets=None,
            initial_state_mean=None, initial_state_covariance=None,
            random_state=None,
            em_vars=['transition_covariance', 'observation_covariance',
                     'initial_state_mean', 'initial_state_covariance'],
            n_dim_state=None, n_dim_obs=None):
        """Initialize Kalman Filter"""
        self.transition_matrices = transition_matrices
        self.observation_matrices = observation_matrices
        self.transition_covariance = transition_covariance
        self.observation_covariance = observation_covariance
        self.transition_offsets = transition_offsets
        self.observation_offsets = observation_offsets
        self.initial_state_mean = initial_state_mean
        self.initial_state_covariance = initial_state_covariance
        self.random_state = random_state
        self.em_vars = em_vars
        self.n_dim_state = n_dim_state
        self.n_dim_obs = n_dim_obs

    def sample(self, n_timesteps, initial_state=None, random_state=None):
        """Sample a state sequence `n_timesteps` timesteps in length.

        Parameters
        ----------
        n_timesteps : int
            number of timesteps

        Returns
        -------
        states : [n_timesteps, n_dim_state] array
            hidden states corresponding to times [0...n_timesteps-1]
        observations : [n_timesteps, n_dim_obs] array
            observations corresponding to times [0...n_timesteps-1]
        """
        (transition_matrices, transition_offsets, transition_covariance,
         observation_matrices, observation_offsets, observation_covariance,
         initial_state_mean, initial_state_covariance) = (
            self._initialize_parameters()
        )

        n_dim_state = transition_matrices.shape[-2]
        n_dim_obs = observation_matrices.shape[-2]
        states = np.zeros((n_timesteps, n_dim_state))
        observations = np.zeros((n_timesteps, n_dim_obs))

        # logic for instantiating rng
        if random_state is None:
            rng = check_random_state(self.random_state)
        else:
            rng = check_random_state(random_state)

        # logic for selecting initial state
        if initial_state is None:
            initial_state = rng.multivariate_normal(
                initial_state_mean,
                initial_state_covariance
            )

        # logic for generating samples
        for t in range(n_timesteps):
            if t == 0:
                states[t] = initial_state
            else:
                transition_matrix = _last_dims(
                    transition_matrices, t - 1
                )
                transition_offset = _last_dims(
                    transition_offsets, t - 1, ndims=1
                )
                transition_covariance = _last_dims(
                    transition_covariance, t - 1
                )
                states[t] = (
                    np.dot(transition_matrix, states[t - 1])
                    + transition_offset
                    + rng.multivariate_normal(
                        np.zeros(n_dim_state),
                        transition_covariance.newbyteorder('=')
                    )
                )

            observation_matrix = _last_dims(
                observation_matrices, t
            )
            observation_offset = _last_dims(
                observation_offsets, t, ndims=1
            )
            observation_covariance = _last_dims(
                observation_covariance, t
            )
            observations[t] = (
                np.dot(observation_matrix, states[t])
                + observation_offset
                + rng.multivariate_normal(
                    np.zeros(n_dim_obs),
                    observation_covariance.newbyteorder('=')
                )
            )

        return (states, np.ma.array(observations))

    def filter(self, X):
        """Apply the Kalman Filter

        Apply the Kalman Filter to estimate the hidden state at time :math:`t`
        for :math:`t = [0...n_timesteps-1]` given observations up to and
        including time `t`.  Observations are assumed to correspond to times
        `[0...n_timesteps-1]`.  The output of this method corresponding to time
        :math:`n_timesteps-1` can be used in :func:`KalmanFilter.filter_update`
        for online updating.

        Parameters
        ----------
        X : [n_timesteps, n_dim_obs] array-like
            observations corresponding to times [0...n_timesteps-1].  If `X` is
            a masked array and any of `X[t]` is masked, then `X[t]` will be
            treated as a missing observation.

        Returns
        -------
        filtered_state_means : [n_timesteps, n_dim_state]
            mean of hidden state distributions for times [0...n_timesteps-1]
            given observations up to and including the current time step
        filtered_state_covariances : [n_timesteps, n_dim_state, n_dim_state] \
        array
            covariance matrix of hidden state distributions for times
            [0...n_timesteps-1] given observations up to and including the
            current time step
        """
        Z = self._parse_observations(X)

        (transition_matrices, transition_offsets, transition_covariance,
         observation_matrices, observation_offsets, observation_covariance,
         initial_state_mean, initial_state_covariance) = (
            self._initialize_parameters()
        )

        (_, _, _, filtered_state_means,
         filtered_state_covariances) = (
            _filter(
                transition_matrices, observation_matrices,
                transition_covariance, observation_covariance,
                transition_offsets, observation_offsets,
                initial_state_mean, initial_state_covariance,
                Z
            )
        )
        return (filtered_state_means, filtered_state_covariances)

    def filter_update(self, filtered_state_mean, filtered_state_covariance,
                      observation=None, transition_matrix=None,
                      transition_offset=None, transition_covariance=None,
                      observation_matrix=None, observation_offset=None,
                      observation_covariance=None):
        r"""Update a Kalman Filter state estimate

        Perform a one-step update to estimate the state at time :math:`t+1`
        give an observation at time :math:`t+1` and the previous estimate for
        time :math:`t` given observations from times :math:`[0...t]`.  This
        method is useful if one wants to track an object with streaming
        observations.

        Parameters
        ----------
        filtered_state_mean : [n_dim_state] array
            mean estimate for state at time t given observations from times
            [1...t]
        filtered_state_covariance : [n_dim_state, n_dim_state] array
            covariance of estimate for state at time t given observations from
            times [1...t]
        observation : [n_dim_obs] array or None
            observation from time t+1.  If `observation` is a masked array and
            any of `observation`'s components are masked or if `observation` is
            None, then `observation` will be treated as a missing observation.
        transition_matrix : optional, [n_dim_state, n_dim_state] array
            state transition matrix from time t to t+1.  If unspecified,
            `self.transition_matrices` will be used.
        transition_offset : optional, [n_dim_state] array
            state offset for transition from time t to t+1.  If unspecified,
            `self.transition_offset` will be used.
        transition_covariance : optional, [n_dim_state, n_dim_state] array
            state transition covariance from time t to t+1.  If unspecified,
            `self.transition_covariance` will be used.
        observation_matrix : optional, [n_dim_obs, n_dim_state] array
            observation matrix at time t+1.  If unspecified,
            `self.observation_matrices` will be used.
        observation_offset : optional, [n_dim_obs] array
            observation offset at time t+1.  If unspecified,
            `self.observation_offset` will be used.
        observation_covariance : optional, [n_dim_obs, n_dim_obs] array
            observation covariance at time t+1.  If unspecified,
            `self.observation_covariance` will be used.

        Returns
        -------
        next_filtered_state_mean : [n_dim_state] array
            mean estimate for state at time t+1 given observations from times
            [1...t+1]
        next_filtered_state_covariance : [n_dim_state, n_dim_state] array
            covariance of estimate for state at time t+1 given observations
            from times [1...t+1]
        """
        # initialize matrices
        (transition_matrices, transition_offsets, transition_covariance,
         observation_matrices, observation_offsets, observation_covariance,
         initial_state_mean, initial_state_covariance) = (
            self._initialize_parameters()
        )
        transition_offset = _arg_or_default(
            transition_offset, self.transition_offsets,
            1, "transition_offset"
        )
        observation_offset = _arg_or_default(
            observation_offset, self.observation_offsets,
            1, "observation_offset"
        )
        transition_matrix = _arg_or_default(
            transition_matrix, self.transition_matrices,
            2, "transition_matrix"
        )
        observation_matrix = _arg_or_default(
            observation_matrix, self.observation_matrices,
            2, "observation_matrix"
        )
        transition_covariance = _arg_or_default(
            transition_covariance, self.transition_covariance,
            2, "transition_covariance"
        )
        observation_covariance = _arg_or_default(
            observation_covariance, self.observation_covariance,
            2, "observation_covariance"
        )

        # Make a masked observation if necessary
        if observation is None:
            n_dim_obs = observation_covariance.shape[0]
            observation = np.ma.array(np.zeros(n_dim_obs))
            observation.mask = True
        else:
            observation = np.ma.asarray(observation)

        predicted_state_mean, predicted_state_covariance = (
            _filter_predict(
                transition_matrix, transition_covariance,
                transition_offset, filtered_state_mean,
                filtered_state_covariance
            )
        )
        (_, next_filtered_state_mean,
         next_filtered_state_covariance) = (
            _filter_correct(
                observation_matrix, observation_covariance,
                observation_offset, predicted_state_mean,
                predicted_state_covariance, observation
            )
        )

        return (next_filtered_state_mean, next_filtered_state_covariance)

    def predict(self, X):
        r"""Apply the Kalman Smoother

        Apply the Kalman Smoother to estimate the hidden state at time
        :math:`t` for :math:`t = [0...n_timesteps-1]` given all observations.
        See :func:`sklearn.kalman._smooth` for more complex output

        Parameters
        ----------
        X : [n_timesteps, n_dim_obs] array-like
            observations corresponding to times [0...n_timesteps-1].  If `X` is
            a masked array and any of `X[t]` is masked, then `X[t]` will be
            treated as a missing observation.

        Returns
        -------
        smoothed_state_means : [n_timesteps, n_dim_state]
            mean of hidden state distributions for times [0...n_timesteps-1]
            given all observations
        """
        Z = self._parse_observations(X)

        (transition_matrices, transition_offsets, transition_covariance,
         observation_matrices, observation_offsets, observation_covariance,
         initial_state_mean, initial_state_covariance) = (
            self._initialize_parameters()
        )

        (predicted_state_means, predicted_state_covariances,
         _, filtered_state_means, filtered_state_covariances) = (
            _filter(
                transition_matrices, observation_matrices,
                transition_covariance, observation_covariance,
                transition_offsets, observation_offsets,
                initial_state_mean, initial_state_covariance, Z
            )
        )
        smoothed_state_means = (
            _smooth(
                transition_matrices, filtered_state_means,
                filtered_state_covariances, predicted_state_means,
                predicted_state_covariances
            )[0]
        )
        return smoothed_state_means

    def fit(self, X, y=None, n_iter=10, em_vars=None):
        """Apply the EM algorithm

        Apply the EM algorithm to estimate all parameters specified by
        `em_vars`.  Note that all variables estimated are assumed to be
        constant for all time.  See :func:`sklearn.kalman._em` for details.

        Parameters
        ----------
        X : [n_timesteps, n_dim_obs] array-like
            observations corresponding to times [0...n_timesteps-1].  If `X` is
            a masked array and any of `X[t]`'s components is masked, then
            `X[t]` will be treated as a missing observation.
        n_iter : int, optional
            number of EM iterations to perform
        em_vars : iterable of strings or 'all'
            variables to perform EM over.  Any variable not appearing here is
            left untouched.
        """
        Z = self._parse_observations(X)

        # initialize parameters
        (self.transition_matrices, self.transition_offsets,
         self.transition_covariance, self.observation_matrices,
         self.observation_offsets, self.observation_covariance,
         self.initial_state_mean, self.initial_state_covariance) = (
            self._initialize_parameters()
        )

        # Create dictionary of variables not to perform EM on
        if em_vars is None:
            em_vars = self.em_vars

        if em_vars == 'all':
            given = {}
        else:
            given = {
                'transition_matrices': self.transition_matrices,
                'observation_matrices': self.observation_matrices,
                'transition_offsets': self.transition_offsets,
                'observation_offsets': self.observation_offsets,
                'transition_covariance': self.transition_covariance,
                'observation_covariance': self.observation_covariance,
                'initial_state_mean': self.initial_state_mean,
                'initial_state_covariance': self.initial_state_covariance
            }
            em_vars = set(em_vars)
            for k in given.keys():
                if k in em_vars:
                    given.pop(k)

        # If a parameter is time varying, print a warning
        for (k, v) in self.get_params().items():
            if k in DIM and (not k in given) and len(v.shape) != DIM[k]:
                warn_str = (
                    '%s has %s dimensions now; after fitting, '
                    + 'it will have dimension %d') % (k, len(v.shape), DIM[k])
                warnings.warn(warn_str)

        # Actual EM iterations
        for i in range(n_iter):
            (predicted_state_means, predicted_state_covariances,
             kalman_gains, filtered_state_means,
             filtered_state_covariances) = (
                _filter(
                    self.transition_matrices, self.observation_matrices,
                    self.transition_covariance, self.observation_covariance,
                    self.transition_offsets, self.observation_offsets,
                    self.initial_state_mean, self.initial_state_covariance,
                    Z
                )
            )
            (smoothed_state_means, smoothed_state_covariances,
             kalman_smoothing_gains) = (
                _smooth(
                    self.transition_matrices, filtered_state_means,
                    filtered_state_covariances, predicted_state_means,
                    predicted_state_covariances
                )
            )
            sigma_pair_smooth = _smooth_pair(
                smoothed_state_covariances,
                kalman_smoothing_gains
            )
            (self.transition_matrices,  self.observation_matrices,
             self.transition_offsets, self.observation_offsets,
             self.transition_covariance, self.observation_covariance,
             self.initial_state_mean, self.initial_state_covariance) = (
                _em(Z, self.transition_offsets, self.observation_offsets,
                    smoothed_state_means, smoothed_state_covariances,
                    sigma_pair_smooth, given=given
                )
            )
        return self

    def score(self, X):
        """Calculate the log likelihood of all observations

        Parameters
        ----------
        X : [n_timesteps, n_dim_obs] array
            observations for time steps [0...n_timesteps-1]

        Returns
        -------
        likelihood : float
            likelihood of all observations
        """
        Z = self._parse_observations(X)

        # initialize parameters
        (transition_matrices, transition_offsets,
         transition_covariance, observation_matrices,
         observation_offsets, observation_covariance,
         initial_state_mean, initial_state_covariance) = (
            self._initialize_parameters()
        )

        # apply the Kalman Filter
        (predicted_state_means, predicted_state_covariances,
         kalman_gains, filtered_state_means,
         filtered_state_covariances) = (
            _filter(
                transition_matrices, observation_matrices,
                transition_covariance, observation_covariance,
                transition_offsets, observation_offsets,
                initial_state_mean, initial_state_covariance,
                Z
            )
        )

        # get likelihoods for each time step
        loglikelihoods = _loglikelihoods(
          observation_matrices, observation_offsets, observation_covariance,
          predicted_state_means, predicted_state_covariances, Z
        )

        return np.sum(loglikelihoods)

    def _initialize_parameters(self):
        """Retrieve parameters if they exist, else replace with defaults"""

        # determine size of state and observation space
        n_dim_state = _determine_dimensionality(
            [(self.transition_matrices, array2d, -2),
             (self.transition_offsets, array1d, -1),
             (self.transition_covariance, array2d, -2),
             (self.initial_state_mean, array1d, -1),
             (self.initial_state_covariance, array2d, -2),
             (self.observation_matrices, array2d, -1)],
            self.n_dim_state
        )
        n_dim_obs = _determine_dimensionality(
            [(self.observation_matrices, array2d, -2),
             (self.observation_offsets, array1d, -1),
             (self.observation_covariance, array2d, -2)],
            self.n_dim_obs
        )

        # initialize parameters
        transition_matrices = (
            np.eye(n_dim_state)
            if self.transition_matrices is None
            else array2d(self.transition_matrices)
        )
        transition_offsets = (
            np.zeros(n_dim_state)
            if self.transition_offsets is None
            else array1d(self.transition_offsets)
        )
        transition_covariance = (
            np.eye(n_dim_state)
            if self.transition_covariance is None
            else array2d(self.transition_covariance)
        )
        observation_matrices = (
            np.eye(n_dim_obs, n_dim_state)
            if self.observation_matrices is None
            else array2d(self.observation_matrices)
        )
        observation_offsets = (
            np.zeros(n_dim_obs)
            if self.observation_offsets is None
            else array1d(self.observation_offsets)
        )
        observation_covariance = (
            np.eye(n_dim_obs)
            if self.observation_covariance is None
            else array2d(self.observation_covariance)
        )
        initial_state_mean = (
            np.zeros(n_dim_state)
            if self.initial_state_mean is None
            else array1d(self.initial_state_mean)
        )
        initial_state_covariance = (
            np.eye(n_dim_state)
            if self.initial_state_covariance is None
            else array2d(self.initial_state_covariance)
        )
        return (transition_matrices, transition_offsets,
                transition_covariance, observation_matrices,
                observation_offsets, observation_covariance,
                initial_state_mean, initial_state_covariance)

    def _parse_observations(self, obs):
        """Safely convert observations to their expected format"""
        obs = np.ma.atleast_2d(obs)
        if obs.shape[0] == 1 and obs.shape[1] > 1:
            obs = obs.T
        return obs
