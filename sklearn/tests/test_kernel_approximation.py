import numpy as np
from scipy.sparse import csr_matrix

from sklearn.utils.testing import assert_array_equal, assert_equal
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_array_almost_equal, assert_raises

from sklearn.metrics.pairwise import kernel_metrics
from sklearn.kernel_approximation import RBFSampler
from sklearn.kernel_approximation import AdditiveChi2Sampler
from sklearn.kernel_approximation import SkewedChi2Sampler
from sklearn.kernel_approximation import Nystroem
from sklearn.metrics.pairwise import polynomial_kernel, rbf_kernel

# generate data
rng = np.random.RandomState(0)
X = rng.random_sample(size=(300, 50))
Y = rng.random_sample(size=(300, 50))
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

    # approximate kernel mapping
    transform = AdditiveChi2Sampler(sample_steps=3)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y)

    kernel_approx = np.dot(X_trans, Y_trans.T)

    assert_array_almost_equal(kernel, kernel_approx, 1)

    X_sp_trans = transform.fit_transform(csr_matrix(X))
    Y_sp_trans = transform.transform(csr_matrix(Y))

    assert_array_equal(X_trans, X_sp_trans.A)
    assert_array_equal(Y_trans, Y_sp_trans.A)

    # test error is raised on negative input
    Y_neg = Y.copy()
    Y_neg[0, 0] = -1
    assert_raises(ValueError, transform.transform, Y_neg)

    # test error on invalid sample_steps
    transform = AdditiveChi2Sampler(sample_steps=4)
    assert_raises(ValueError, transform.fit, X)

    # test that the sample interval is set correctly
    sample_steps_available = [1, 2, 3]
    for sample_steps in sample_steps_available:

        # test that the sample_interval is initialized correctly
        transform = AdditiveChi2Sampler(sample_steps=sample_steps)
        assert_equal(transform.sample_interval, None)

        # test that the sample_interval is changed in the fit method
        transform.fit(X)
        assert_not_equal(transform.sample_interval_, None)

    # test that the sample_interval is set correctly
    sample_interval = 0.3
    transform = AdditiveChi2Sampler(sample_steps=4,
                                    sample_interval=sample_interval)
    assert_equal(transform.sample_interval, sample_interval)
    transform.fit(X)
    assert_equal(transform.sample_interval_, sample_interval)


def test_skewed_chi2_sampler():
    """test that RBFSampler approximates kernel on random data"""

    # compute exact kernel
    c = 0.03
    # appreviations for easier formular
    X_c = (X + c)[:, np.newaxis, :]
    Y_c = (Y + c)[np.newaxis, :, :]

    # we do it in log-space in the hope that it's more stable
    # this array is n_samples_x x n_samples_y big x n_features
    log_kernel = ((np.log(X_c) / 2.) + (np.log(Y_c) / 2.) + np.log(2.) -
                  np.log(X_c + Y_c))
    # reduce to n_samples_x x n_samples_y by summing over features in log-space
    kernel = np.exp(log_kernel.sum(axis=2))

    # approximate kernel mapping
    transform = SkewedChi2Sampler(skewedness=c, n_components=1000,
                                  random_state=42)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y)

    kernel_approx = np.dot(X_trans, Y_trans.T)
    assert_array_almost_equal(kernel, kernel_approx, 1)

    # test error is raised on negative input
    Y_neg = Y.copy()
    Y_neg[0, 0] = -1
    assert_raises(ValueError, transform.transform, Y_neg)


def test_rbf_sampler():
    """test that RBFSampler approximates kernel on random data"""
    # compute exact kernel
    gamma = 10.
    kernel = rbf_kernel(X, Y, gamma=gamma)

    # approximate kernel mapping
    rbf_transform = RBFSampler(gamma=gamma, n_components=1000, random_state=42)
    X_trans = rbf_transform.fit_transform(X)
    Y_trans = rbf_transform.transform(Y)
    kernel_approx = np.dot(X_trans, Y_trans.T)

    assert_array_almost_equal(kernel, kernel_approx, 1)


def test_input_validation():
    """Regression test: kernel approx. transformers should work on lists

    No assertions; the old versions would simply crash
    """
    X = [[1, 2], [3, 4], [5, 6]]
    AdditiveChi2Sampler().fit(X).transform(X)
    SkewedChi2Sampler().fit(X).transform(X)
    RBFSampler().fit(X).transform(X)

    X = csr_matrix(X)
    RBFSampler().fit(X).transform(X)


def test_nystroem_approximation_with_number_samples_is_exact():
    # some basic tests
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 4))

    # With n_components = n_samples this is exact
    ny_random = Nystroem(n_components=X.shape[0], basis_method='random')
    X_transformed_random = ny_random.fit_transform(X)
    K = rbf_kernel(X)
    assert_array_equal(np.sort(ny_random.component_indices_), np.arange(X.shape[0]))
    assert_array_almost_equal(np.dot(X_transformed_random, X_transformed_random.T), K)

    ny_clustered = Nystroem(n_components=X.shape[0], basis_method='clustered')
    X_transformed_clustered = ny_clustered.fit_transform(X)
    K = rbf_kernel(X)
    # No component indicies to report for k-means
    assert_equal(ny_clustered.component_indices_, None)
    assert_array_almost_equal(np.dot(X_transformed_clustered, X_transformed_clustered.T), K)


def test_nystroem_approximation_returns_appropriate_indices():
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 4))

    ny_random = Nystroem(n_components=2, basis_method='random')
    X_transformed = ny_random.fit_transform(X)
    assert_equal(X_transformed.shape, (X.shape[0], 2))
    assert_equal(len(ny_random.component_indices_), 2)
    assert_array_almost_equal(ny_random.components_, X[ny_random.component_indices_])

    ny_clustered = Nystroem(n_components=2, basis_method='clustered')
    ny_clustered.fit_transform(X)
    # No component indicies to report for k-means
    assert_equal(ny_clustered.component_indices_, None)


def test_nystroem_approximation_with_singular_kernel_matrix():
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 4))
    X = np.concatenate((X, X[-2:, :]), axis=0)

    K = rbf_kernel(X)
    assert_equal(np.linalg.matrix_rank(K), 10)

    ny_random = Nystroem(n_components=X.shape[0], basis_method='random')
    X_transformed = ny_random.fit_transform(X)
    assert_equal(X_transformed.shape, (X.shape[0], 12))
    assert_array_almost_equal(np.dot(X_transformed, X_transformed.T), K)


def test_nystroem_approximation_for_multiple_kernels():
    """test that Nystroem approximates kernel on random data"""
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 4))
    trans_not_valid = Nystroem(n_components=2, random_state=rnd,
                               basis_method="not_a_valid_basis_method")
    assert_raises(NameError, trans_not_valid.fit, X)

    # Kernel tests to perform with each basis method used
    def test_nystroem_approximation_with_basis(tested_basis):
         # Test default kernel
        trans = Nystroem(n_components=2, random_state=rnd, basis_method=tested_basis)
        transformed = trans.fit(X).transform(X)
        assert_equal(transformed.shape, (X.shape[0], 2))

        # test callable kernel
        linear_kernel = lambda X, Y: np.dot(X, Y.T)
        trans = Nystroem(n_components=2, kernel=linear_kernel, random_state=rnd, basis_method=tested_basis)
        transformed = trans.fit(X).transform(X)
        assert_equal(transformed.shape, (X.shape[0], 2))

        # test that available kernels fit and transform
        kernels_available = kernel_metrics()
        for kern in kernels_available:
            trans = Nystroem(n_components=2, kernel=kern, random_state=rnd, basis_method=tested_basis)
            transformed = trans.fit(X).transform(X)
            assert_equal(transformed.shape, (X.shape[0], 2))

        # Test default kernel
        trans = Nystroem(n_components=2, random_state=rnd, basis_method=tested_basis)
        transformed = trans.fit(X).transform(X)
        assert_equal(transformed.shape, (X.shape[0], 2))

        # test callable kernel
        linear_kernel = lambda X, Y: np.dot(X, Y.T)
        trans = Nystroem(n_components=2, kernel=linear_kernel, random_state=rnd, basis_method=tested_basis)
        transformed = trans.fit(X).transform(X)
        assert_equal(transformed.shape, (X.shape[0], 2))

        # test that available kernels fit and transform
        kernels_available = kernel_metrics()
        for kern in kernels_available:
            trans = Nystroem(n_components=2, kernel=kern, random_state=rnd, basis_method=tested_basis)
            transformed = trans.fit(X).transform(X)
            assert_equal(transformed.shape, (X.shape[0], 2))

    # Go through all the kernels with each basis_method
    basis_methods = ("random", "clustered")
    for current_basis in basis_methods:
        yield test_nystroem_approximation_with_basis, current_basis


def test_nystroem_poly_kernel_params():
    """Non-regression: Nystroem should pass other parameters beside gamma."""
    rnd = np.random.RandomState(37)
    X = rnd.uniform(size=(10, 4))

    K = polynomial_kernel(X, degree=3.1, coef0=.1)
    nystroem_random = Nystroem(kernel="polynomial", n_components=X.shape[0],
                               degree=3.1, coef0=.1, basis_method="random")
    nystroem_k_means = Nystroem(kernel="polynomial", n_components=X.shape[0],
                                degree=3.1, coef0=.1, basis_method="clustered")

    transformed_k_means = nystroem_k_means.fit_transform(X)
    transformed_random = nystroem_random.fit_transform(X)

    assert_array_almost_equal(np.dot(transformed_k_means,
                                     transformed_k_means.T), K)
    assert_array_almost_equal(np.dot(transformed_random,
                                     transformed_random.T), K)


def test_nystroem_callable():
    """Test Nystroem on a callable."""
    rnd = np.random.RandomState(42)
    n_samples = 10
    X = rnd.uniform(size=(n_samples, 4))

    def logging_histogram_kernel(x, y, log):
        """Histogram kernel that writes to a log."""
        log.append(1)
        return np.minimum(x, y).sum()

    kernel_log = []
    Nystroem(kernel=logging_histogram_kernel,
             n_components=(n_samples - 1),
             kernel_params={'log': kernel_log}, basis_method="clustered").fit(X)

    assert_equal(len(kernel_log), n_samples * (n_samples - 1) / 2)

    kernel_log = []
    Nystroem(kernel=logging_histogram_kernel,
             n_components=(n_samples - 1),
             kernel_params={'log': kernel_log}, basis_method="random").fit(X)

    assert_equal(len(kernel_log), n_samples * (n_samples - 1) / 2)
