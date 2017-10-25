import numpy as np
from numpy.testing import assert_array_almost_equal

from sklearn.utils.estimator_checks import check_estimator
from ..sfa import SFA


def _input_data_sin_waves(n_samples):
    """ Support function to create two mixes sine waves to test SFA.

    By construction, the first sine wave is faster (higher frequency) than
    the second.

    Sine waves are the analytical solution to the SFA problem
    (Wiskott and Sejnowski, 2002), so we are sure that SFA is going to
    be able to recover them.
    """
    t = np.linspace(0.0, 1.0, num=n_samples)
    X = np.array([np.sin(2 * np.pi * 5 * t),
                  np.sin(2 * np.pi * t)]).T

    X = X - X[:-1, :].mean(axis=0)
    X = X / X[:-1, :].std(axis=0)

    transformation_mix = np.array([[0.3, 0.9], [0.4, 0.1]])
    transformation_mean = np.array([10.3, -29.99])
    X_mixed = np.dot(X, transformation_mix) + transformation_mean

    return X, X_mixed


def test_sfa_is_valid_estimator():
    return check_estimator(SFA)


def test_sfa():
    # Create two sin waves with zero mean and unit variance, and mix them
    # together with a arbitrary linear transform. SFA should be able to
    # unmix them.
    n_samples = 10000
    X, X_mixed = _input_data_sin_waves(n_samples)

    sfa = SFA()
    sfa.fit(X_mixed)
    X_output = sfa.transform(X_mixed)

    # Output should have zero mean, unit variance
    assert_array_almost_equal(X_output.mean(axis=0), 0)
    assert_array_almost_equal(X_output.std(axis=0), 1)

    # The output signals are perfectly correlated with the original,
    # unmixed signals, ordered by slowness (slowest first).
    correlation = np.dot(X[:, [1, 0]].T, X_output) / n_samples
    assert_array_almost_equal(abs(correlation), np.eye(2), 4)


def test_sfa_reduce_output_dim():
    n_samples = 10000
    X, X_mixed = _input_data_sin_waves(n_samples)

    sfa = SFA(n_components=1)
    sfa.fit(X_mixed)
    X_output = sfa.transform(X_mixed)

    assert X_output.shape[1] == 1

    # The output signals is *the slowest* sine wave.
    correlation = np.dot(X[:, 1:].T, X_output) / n_samples
    assert_array_almost_equal(abs(correlation), 1.0, 4)
