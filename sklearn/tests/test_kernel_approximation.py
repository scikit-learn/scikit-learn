import numpy as np
from scipy.spatial.distance import cdist

from ..kernel_approximation import RBFSampler
from ..kernel_approximation import AdditiveChi2Sampler
from ..kernel_approximation import SkewedChi2Sampler

# generate data
X = np.random.uniform(size=(300, 50))
Y = np.random.uniform(size=(300, 50))
X /= X.sum(axis=1)[:, np.newaxis]
Y /= Y.sum(axis=1)[:, np.newaxis]


def test_additive_chi2_sampler():
    """test that AdditiveChi2Sampler approximates kernel on random data"""

    # compute exact kernel
    # appreviations for easier formular
    X_ = X[:, np.newaxis, :]
    Y_ = Y[np.newaxis, :, :]

    large_kernel = 2 * X_ * Y_ / (X_ + Y_)

    # reduce to n_samples_x x n_samples_y by summing over features
    kernel = (large_kernel.sum(axis=2))

    # appoximate kernel mapping
    transform = AdditiveChi2Sampler(sample_steps=3)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y)
    kernel_approx = np.dot(X_trans, Y_trans.T)

    np.testing.assert_array_almost_equal(kernel, kernel_approx, 1)


def test_skewed_chi2_sampler():
    """test that RBFSampler approximates kernel on random data"""

    # compute exact kernel
    c = 0.03
    # appreviations for easier formular
    X_c = (X + c)[:, np.newaxis, :]
    Y_c = (Y + c)[np.newaxis, :, :]

    # we do it in log-space in the hope that it's more stable
    # this array is n_samples_x x n_samples_y big x n_features
    log_kernel = ((np.log(X_c) / 2.) + (np.log(Y_c) / 2.) +
        np.log(2.) - np.log(X_c + Y_c))
    # reduce to n_samples_x x n_samples_y by summing over features in log-space
    kernel = np.exp(log_kernel.sum(axis=2))

    # appoximate kernel mapping
    transform = SkewedChi2Sampler(skewedness=c, n_components=1000, random_state=42)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y)
    kernel_approx = np.dot(X_trans, Y_trans.T)
    np.testing.assert_array_almost_equal(kernel, kernel_approx, 1)

def test_rbf_sampler():
    """test that RBFSampler approximates kernel on random data"""
    # compute exact kernel
    gamma = 10.
    dists = cdist(X, Y)
    kernel = np.exp(-gamma * dists ** 2)

    # appoximate kernel mapping
    rbf_transform = RBFSampler(gamma=gamma, n_components=1000, random_state=42)
    X_trans = rbf_transform.fit_transform(X)
    Y_trans = rbf_transform.transform(Y)
    kernel_approx = np.dot(X_trans, Y_trans.T)

    np.testing.assert_array_almost_equal(kernel, kernel_approx, 1)

if __name__ == "__main__":
    test_additive_chi2_sampler()
    test_skewed_chi2_sampler()
    test_rbf_sampler()
