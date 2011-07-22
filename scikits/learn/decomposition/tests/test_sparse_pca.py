# Author: Vlad Niculae
# License: BSD

import sys

import numpy as np
from .. import SparsePCA, MiniBatchSparsePCA, dict_learning_online
from ..sparse_pca import _update_code, _update_code_parallel


from numpy.testing import assert_array_almost_equal, assert_equal


def generate_toy_data(n_atoms, n_samples, image_size):
    n_features = image_size[0] * image_size[1]

    np.random.seed(0)
    U = np.random.randn(n_samples, n_atoms)
    V = np.random.randn(n_atoms, n_features)

    centers = [(3, 3), (6, 7), (8, 1)]
    sz = [1, 2, 1]
    for k in range(n_atoms):
        img = np.zeros(image_size)
        xmin, xmax = centers[k][0] - sz[k], centers[k][0] + sz[k]
        ymin, ymax = centers[k][1] - sz[k], centers[k][1] + sz[k]
        img[xmin:xmax][:, ymin:ymax] = 1.0
        V[k, :] = img.ravel()

    # Y is defined by : Y = UV + noise
    Y = np.dot(U, V)
    Y += 0.1 * np.random.randn(Y.shape[0], Y.shape[1])  # Add noise
    return Y, U, V

# SparsePCA can be a bit slow. To avoid having test times go up, we
# test different aspects of the code in the same test


def test_correct_shapes():
    np.random.seed(0)
    X = np.random.randn(12, 10)
    pca = SparsePCA(n_components=8)
    U = pca.fit_transform(X)
    assert_equal(pca.components_.shape, (8, 10))
    assert_equal(U.shape, (12, 8))
    # test overcomplete decomposition
    pca = SparsePCA(n_components=13)
    U = pca.fit_transform(X)
    assert_equal(pca.components_.shape, (13, 10))
    assert_equal(U.shape, (12, 13))


def test_fit_transform():
    Y, _, _ = generate_toy_data(3, 10, (8, 8))  # wide array
    spca_lars = SparsePCA(n_components=3, method='lars').fit(Y)
    U1 = spca_lars.transform(Y)
    # Test multiple CPUs
    if sys.platform == 'win32':  # fake parallelism for win32
        import scikits.learn.externals.joblib.parallel as joblib_par
        _mp = joblib_par.multiprocessing
        joblib_par.multiprocessing = None
        try:
            U2 = SparsePCA(n_components=3, n_jobs=2).fit(Y).transform(Y)
        finally:
            joblib_par.multiprocessing = _mp
    else:  # we can efficiently use parallelism
        U2 = SparsePCA(n_components=3, n_jobs=2).fit(Y).transform(Y)
    assert_array_almost_equal(U1, U2)
    # Test that CD gives similar results
    spca_lasso = SparsePCA(n_components=3, method='cd').fit(Y)
    assert_array_almost_equal(spca_lasso.components_, spca_lars.components_)


def test_fit_transform_tall():
    Y, _, _ = generate_toy_data(3, 65, (8, 8))  # tall array
    U1 = SparsePCA(n_components=3, method='lars').fit_transform(Y)
    U2 = SparsePCA(n_components=3, method='cd').fit(Y).transform(Y)
    assert_array_almost_equal(U1, U2)


def test_sparse_code():
    np.random.seed(0)
    dictionary = np.random.randn(10, 3)
    real_code = np.zeros((3, 5))
    real_code.ravel()[np.random.randint(15, size=6)] = 1.0
    Y = np.dot(dictionary, real_code)
    est_code_1 = _update_code(dictionary, Y, alpha=1.0)
    est_code_2 = _update_code_parallel(dictionary, Y, alpha=1.0)
    assert_equal(est_code_1.shape, real_code.shape)
    assert_equal(est_code_1, est_code_2)
    assert_equal(est_code_1.nonzero(), real_code.nonzero())


def test_initialization():
    U_init = np.random.randn(5, 3)
    V_init = np.random.randn(3, 4)
    model = SparsePCA(n_components=3, U_init=U_init, V_init=V_init, max_iter=0)
    model.fit(np.random.randn(5, 4))
    assert_equal(model.components_, V_init)


def test_dict_learning_online_shapes():
    np.random.seed(0)
    X = np.random.randn(12, 10)
    codeT, dictionaryT = dict_learning_online(X.T, n_atoms=8, alpha=1)
    assert_equal(codeT.shape, (8, 12))
    assert_equal(dictionaryT.shape, (10, 8))
    assert_equal(np.dot(codeT.T, dictionaryT.T).shape, X.shape)


def test_mini_batch_sparse_pca_correct_shapes():
    np.random.seed(0)
    X = np.random.randn(12, 10)
    pca = MiniBatchSparsePCA(n_components=8)
    U = pca.fit_transform(X)
    assert_equal(pca.components_.shape, (8, 10))
    assert_equal(U.shape, (12, 8))
    # test overcomplete decomposition
    pca = MiniBatchSparsePCA(n_components=13)
    U = pca.fit_transform(X)
    assert_equal(pca.components_.shape, (13, 10))
    assert_equal(U.shape, (12, 13))