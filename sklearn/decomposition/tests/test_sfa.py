import numpy as np
from numpy.testing import assert_array_almost_equal

from sklearn.utils.estimator_checks import check_estimator
from ..sfa import SFA


def test_sfa_is_valid_estimator():
    return check_estimator(SFA)


def test_sfa():
    # Create two sin waves with zero mean and unit variance, and mix them
    # together with a arbitrary linear transform. SFA should be able to
    # unmix them.
    n_samples = 10000
    t = np.linspace(0.0, 1.0, num=n_samples)
    X = np.array([np.sin(2 * np.pi * t),
                  np.sin(2 * np.pi * 5 * t)]).T
    X = X - X[:-1, :].mean(axis=0)
    X = X / X[:-1, :].std(axis=0)

    transformation_mix = np.array([[0.3, 0.9], [0.4, 0.1]])
    transformation_mean = np.array([10.3, -29.99])
    X_mixed = np.dot(X, transformation_mix) + transformation_mean

    sfa = SFA()
    sfa.fit(X_mixed)
    X_output = sfa.transform(X_mixed)

    # Output should have zero mean, unit variance
    assert_array_almost_equal(X_output.mean(axis=0), 0)
    assert_array_almost_equal(X_output.std(axis=0), 1)

    # The output signals are perfectly correlated with the original,
    # unmixed signals.
    correlation = np.dot(X.T, X_output) / n_samples
    assert_array_almost_equal(abs(correlation), np.eye(2), 4)
