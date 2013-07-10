"""Testing for Spectral Biclustering methods"""

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import SkipTest

from sklearn.bicluster.spectral import SpectralCoclustering
from sklearn.bicluster.spectral import SpectralBiclustering
from sklearn.bicluster.spectral import _sparse_min
from sklearn.bicluster.spectral import _scale_preprocess
from sklearn.bicluster.spectral import _bistochastic_preprocess
from sklearn.bicluster.spectral import _log_preprocess

from sklearn.datasets import make_biclusters, make_checkerboard

import numpy as np
from scipy.sparse import csr_matrix, issparse
from itertools import permutations


def _check_label_permutations(a, b, n_labels):
    # TODO: replace with hungarian algorithm, when it is implemented
    assert(np.any(np.array(p)[a] == b
                  for p in permutations(range(n_labels))))


def test_spectral_coclustering():
    """Test Dhillon's Spectral CoClustering on a simple problem."""
    random_state = 0
    for noise in (0, 0.5):
        S, rows, cols = make_biclusters((30, 30), 3, noise=noise,
                                        random_state=random_state)
        for mat in (S, csr_matrix(S)):
            for svd_method in ('randomized', 'arpack'):
                model = SpectralCoclustering(n_clusters=3,
                                             svd_method=svd_method,
                                             random_state=random_state)
                model.fit(mat)

                assert_equal(model.rows_.shape, (3, 30))
                assert_array_equal(model.rows_.sum(axis=0), np.ones(30))
                assert_array_equal(model.columns_.sum(axis=0), np.ones(30))
                _check_label_permutations(model.rows_, rows, 3)
                _check_label_permutations(model.columns_, cols, 3)


def test_spectral_biclustering():
    """Test Kluger methods on a checkerboard dataset."""
    raise SkipTest('Permutations are slow. Skipping until'
                   ' faster matching algorithm is available.')
    random_state = 0
    for noise in (0.5, 0):
        S, rows, cols = make_checkerboard((30, 30), 3, noise=noise,
                                          random_state=random_state)
        for n_clusters in ((3, 3), 3):
            for mat in (S, csr_matrix(S)):
                for method in ('scale', 'bistochastic', 'log'):
                    for svd_method in ('randomized', 'arpack'):
                        if svd_method == 'arpack':
                            # fails with default value
                            svd_kwargs = {'ncv' : 20}
                        else:
                            svd_kwargs = {}
                        model = SpectralBiclustering(n_clusters=n_clusters,
                                                     method=method,
                                                     svd_method=svd_method,
                                                     svd_kwargs=svd_kwargs,
                                                     random_state=random_state)
                        model.fit(mat)

                        assert_equal(model.rows_.shape, (9, 30))
                        assert_equal(model.columns_.shape, (9, 30))
                        assert_array_equal(model.rows_.sum(axis=0), np.repeat(3, 30))
                        assert_array_equal(model.columns_.sum(axis=0), np.repeat(3, 30))
                        _check_label_permutations(model.rows_, rows, 3)
                        _check_label_permutations(model.columns_, cols, 3)


def _do_scale_test(scaled):
    """Check that rows sum to one constant, and columns to another."""
    row_sum = scaled.sum(axis=1)
    col_sum = scaled.sum(axis=0)
    if issparse(scaled):
        row_sum = np.asarray(row_sum).squeeze()
        col_sum = np.asarray(col_sum).squeeze()
    assert_array_almost_equal(row_sum, np.tile(row_sum.mean(), 100),
                              decimal=1)
    assert_array_almost_equal(col_sum, np.tile(col_sum.mean(), 100),
                              decimal=1)


def _do_bistochastic_test(scaled):
    """Check that rows and columns sum to the same constant."""
    _do_scale_test(scaled)
    assert_almost_equal(scaled.sum(axis=0).mean(),
                        scaled.sum(axis=1).mean(),
                        decimal=1)


def test_scale_preprocess():
    generator = np.random.RandomState(0)
    X = generator.rand(100, 100)
    for mat in (X, csr_matrix(X)):
        scaled, _, _ = _scale_preprocess(mat)
        _do_scale_test(scaled)
        if issparse(mat):
            assert issparse(scaled)


def test_bistochastic_preprocess():
    generator = np.random.RandomState(0)
    X = generator.rand(100, 100)
    for mat in (X, csr_matrix(X)):
        scaled = _bistochastic_preprocess(mat)
        _do_bistochastic_test(scaled)
        if issparse(mat):
            assert issparse(scaled)


def test_log_preprocess():
    # adding any constant to a log-scaled matrix should make it
    # bistochastic
    generator = np.random.RandomState(0)
    X = generator.rand(100, 100)
    for mat in (X, csr_matrix(X)):
        scaled = _log_preprocess(mat) + 1
        _do_bistochastic_test(scaled)
        assert not issparse(scaled)


def test_fit_best_piecewise():
    model = SpectralBiclustering(random_state=0)
    vectors = np.array([[0, 0, 0, 1, 1, 1],
                        [2, 2, 2, 3, 3, 3],
                        [0, 1, 2, 3, 4, 5],])
    best = model._fit_best_piecewise(vectors, n_best=2, n_clusters=2)
    assert_array_equal(best, vectors[:2])


def test_project_and_cluster():
    model = SpectralBiclustering(random_state=0)
    data = np.array([[1, 1, 1,],
                     [1, 1, 1,],
                     [3, 6, 3,],
                     [3, 6, 3,]])
    vectors = np.array([[1, 0,],
                        [0, 1,],
                        [0, 0,]])
    for mat in (data, csr_matrix(data)):
        labels = model._project_and_cluster(data, vectors,
                                            n_clusters=2)
        assert_array_equal(labels, [0, 0, 1, 1])


def test_errors():
    data = np.arange(25).reshape((5, 5))

    model = SpectralBiclustering(n_clusters=(3, 3, 3))
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering(n_clusters='abc')
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering(n_clusters=(3, 'abc'))
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering(method='unknown')
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering(svd_method='unknown')
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering(n_components=3, n_best=4)
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering()
    data = np.arange(27).reshape((3, 3, 3))
    assert_raises(ValueError, model.fit, data)
