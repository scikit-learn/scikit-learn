# Author: Vlad Niculae
# License: BSD 3 clause

import sys
import pytest

import numpy as np

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import if_safe_multiprocessing_with_blas

from sklearn.decomposition import SparsePCA, MiniBatchSparsePCA, PCA
from sklearn.utils import check_random_state


def generate_toy_data(n_components, n_samples, image_size, random_state=None):
    n_features = image_size[0] * image_size[1]

    rng = check_random_state(random_state)
    U = rng.randn(n_samples, n_components)
    V = rng.randn(n_components, n_features)

    centers = [(3, 3), (6, 7), (8, 1)]
    sz = [1, 2, 1]
    for k in range(n_components):
        img = np.zeros(image_size)
        xmin, xmax = centers[k][0] - sz[k], centers[k][0] + sz[k]
        ymin, ymax = centers[k][1] - sz[k], centers[k][1] + sz[k]
        img[xmin:xmax][:, ymin:ymax] = 1.0
        V[k, :] = img.ravel()

    # Y is defined by : Y = UV + noise
    Y = np.dot(U, V)
    Y += 0.1 * rng.randn(Y.shape[0], Y.shape[1])  # Add noise
    return Y, U, V

# SparsePCA can be a bit slow. To avoid having test times go up, we
# test different aspects of the code in the same test


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
def test_correct_shapes(norm_comp):
    rng = np.random.RandomState(0)
    X = rng.randn(12, 10)
    spca = SparsePCA(n_components=8, random_state=rng,
                     normalize_components=norm_comp)
    U = spca.fit_transform(X)
    assert_equal(spca.components_.shape, (8, 10))
    assert_equal(U.shape, (12, 8))
    # test overcomplete decomposition
    spca = SparsePCA(n_components=13, random_state=rng,
                     normalize_components=norm_comp)
    U = spca.fit_transform(X)
    assert_equal(spca.components_.shape, (13, 10))
    assert_equal(U.shape, (12, 13))


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
def test_fit_transform(norm_comp):
    alpha = 1
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 10, (8, 8), random_state=rng)  # wide array
    spca_lars = SparsePCA(n_components=3, method='lars', alpha=alpha,
                          random_state=0, normalize_components=norm_comp)
    spca_lars.fit(Y)

    # Test that CD gives similar results
    spca_lasso = SparsePCA(n_components=3, method='cd', random_state=0,
                           alpha=alpha, normalize_components=norm_comp)
    spca_lasso.fit(Y)
    assert_array_almost_equal(spca_lasso.components_, spca_lars.components_)

    # Test that deprecated ridge_alpha parameter throws warning
    warning_msg = "The ridge_alpha parameter on transform()"
    assert_warns_message(DeprecationWarning, warning_msg, spca_lars.transform,
                         Y, ridge_alpha=0.01)
    assert_warns_message(DeprecationWarning, warning_msg, spca_lars.transform,
                         Y, ridge_alpha=None)


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
@if_safe_multiprocessing_with_blas
def test_fit_transform_parallel(norm_comp):
    alpha = 1
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 10, (8, 8), random_state=rng)  # wide array
    spca_lars = SparsePCA(n_components=3, method='lars', alpha=alpha,
                          random_state=0, normalize_components=norm_comp)
    spca_lars.fit(Y)
    U1 = spca_lars.transform(Y)
    # Test multiple CPUs
    spca = SparsePCA(n_components=3, n_jobs=2, method='lars', alpha=alpha,
                     random_state=0, normalize_components=norm_comp).fit(Y)
    U2 = spca.transform(Y)
    assert_true(not np.all(spca_lars.components_ == 0))
    assert_array_almost_equal(U1, U2)


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
def test_transform_nan(norm_comp):
    # Test that SparsePCA won't return NaN when there is 0 feature in all
    # samples.
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 10, (8, 8), random_state=rng)  # wide array
    Y[:, 0] = 0
    estimator = SparsePCA(n_components=8, normalize_components=norm_comp)
    assert_false(np.any(np.isnan(estimator.fit_transform(Y))))


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
def test_fit_transform_tall(norm_comp):
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 65, (8, 8), random_state=rng)  # tall array
    spca_lars = SparsePCA(n_components=3, method='lars',
                          random_state=rng, normalize_components=norm_comp)
    U1 = spca_lars.fit_transform(Y)
    spca_lasso = SparsePCA(n_components=3, method='cd',
                           random_state=rng, normalize_components=norm_comp)
    U2 = spca_lasso.fit(Y).transform(Y)
    assert_array_almost_equal(U1, U2)


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
def test_initialization(norm_comp):
    rng = np.random.RandomState(0)
    U_init = rng.randn(5, 3)
    V_init = rng.randn(3, 4)
    model = SparsePCA(n_components=3, U_init=U_init, V_init=V_init, max_iter=0,
                      random_state=rng, normalize_components=norm_comp)
    model.fit(rng.randn(5, 4))
    if norm_comp:
        assert_allclose(model.components_,
                        V_init / np.linalg.norm(V_init, axis=1)[:, None])
    else:
        assert_allclose(model.components_, V_init)


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
def test_mini_batch_correct_shapes(norm_comp):
    rng = np.random.RandomState(0)
    X = rng.randn(12, 10)
    pca = MiniBatchSparsePCA(n_components=8, random_state=rng,
                             normalize_components=norm_comp)
    U = pca.fit_transform(X)
    assert_equal(pca.components_.shape, (8, 10))
    assert_equal(U.shape, (12, 8))
    # test overcomplete decomposition
    pca = MiniBatchSparsePCA(n_components=13, random_state=rng,
                             normalize_components=norm_comp)
    U = pca.fit_transform(X)
    assert_equal(pca.components_.shape, (13, 10))
    assert_equal(U.shape, (12, 13))


@pytest.mark.filterwarnings("ignore:normalize_components")
@pytest.mark.parametrize("norm_comp", [False, True])
def test_mini_batch_fit_transform(norm_comp):
    raise SkipTest("skipping mini_batch_fit_transform.")
    alpha = 1
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 10, (8, 8), random_state=rng)  # wide array
    spca_lars = MiniBatchSparsePCA(n_components=3, random_state=0,
                                   alpha=alpha,
                                   normalize_components=norm_comp).fit(Y)
    U1 = spca_lars.transform(Y)
    # Test multiple CPUs
    if sys.platform == 'win32':  # fake parallelism for win32
        import sklearn.utils._joblib.parallel as joblib_par
        _mp = joblib_par.multiprocessing
        joblib_par.multiprocessing = None
        try:
            spca = MiniBatchSparsePCA(n_components=3, n_jobs=2, alpha=alpha,
                                      random_state=0,
                                      normalize_components=norm_comp)
            U2 = spca.fit(Y).transform(Y)
        finally:
            joblib_par.multiprocessing = _mp
    else:  # we can efficiently use parallelism
        spca = MiniBatchSparsePCA(n_components=3, n_jobs=2, alpha=alpha,
                                  random_state=0,
                                  normalize_components=norm_comp)
        U2 = spca.fit(Y).transform(Y)
    assert_true(not np.all(spca_lars.components_ == 0))
    assert_array_almost_equal(U1, U2)
    # Test that CD gives similar results
    spca_lasso = MiniBatchSparsePCA(n_components=3, method='cd', alpha=alpha,
                                    random_state=0,
                                    normalize_components=norm_comp).fit(Y)
    assert_array_almost_equal(spca_lasso.components_, spca_lars.components_)


def test_scaling_fit_transform():
    alpha = 1
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 1000, (8, 8), random_state=rng)
    spca_lars = SparsePCA(n_components=3, method='lars', alpha=alpha,
                          random_state=rng, normalize_components=True)
    results_train = spca_lars.fit_transform(Y)
    results_test = spca_lars.transform(Y[:10])
    assert_allclose(results_train[0], results_test[0])


def test_pca_vs_spca():
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 1000, (8, 8), random_state=rng)
    Z, _, _ = generate_toy_data(3, 10, (8, 8), random_state=rng)
    spca = SparsePCA(alpha=0, ridge_alpha=0, n_components=2,
                     normalize_components=True)
    pca = PCA(n_components=2)
    pca.fit(Y)
    spca.fit(Y)
    results_test_pca = pca.transform(Z)
    results_test_spca = spca.transform(Z)
    assert_allclose(np.abs(spca.components_.dot(pca.components_.T)),
                    np.eye(2), atol=1e-5)
    results_test_pca *= np.sign(results_test_pca[0, :])
    results_test_spca *= np.sign(results_test_spca[0, :])
    assert_allclose(results_test_pca, results_test_spca)


@pytest.mark.parametrize("spca", [SparsePCA, MiniBatchSparsePCA])
def test_spca_deprecation_warning(spca):
    rng = np.random.RandomState(0)
    Y, _, _ = generate_toy_data(3, 10, (8, 8), random_state=rng)
    warn_message = "normalize_components"
    assert_warns_message(DeprecationWarning, warn_message,
                         spca(normalize_components=False).fit, Y)
