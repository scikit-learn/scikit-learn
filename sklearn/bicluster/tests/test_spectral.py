"""Testing for Spectral Biclustering methods"""

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises

from sklearn.bicluster.spectral import SpectralCoclustering
from sklearn.bicluster.spectral import SpectralBiclustering
from sklearn.bicluster.spectral import _scale_preprocess
from sklearn.bicluster.spectral import _bistochastic_preprocess
from sklearn.bicluster.spectral import _log_preprocess
from sklearn.bicluster.spectral import _fit_best_piecewise
from sklearn.bicluster.spectral import _project_and_cluster

from sklearn.datasets import make_biclusters, make_checkerboard

import numpy as np

from itertools import permutations


def _check_label_permutations(a, b, n_labels):
    # TODO: replace with hungarian algorithm, when it is implemented
    assert(np.any(np.array(p)[a] == b
                  for p in permutations(range(n_labels))))


def test_spectral_biclustering_dhillon():
    """Test Dhillon's Spectral CoClustering on a simple problem."""
    random_state = 0
    for noise in (0, 0.5):
        S, rows, cols = make_biclusters((30, 30), 3, noise=noise,
                                        random_state=random_state)
        for svd_method in ('randomized', 'arpack'):
            model = SpectralCoclustering(n_clusters=3,
                                         svd_method=svd_method,
                                         random_state=random_state)
            model.fit(S)

            assert_equal(model.rows_.shape, (3, 30))
            assert_array_equal(model.rows_.sum(axis=0), np.ones(30))
            assert_array_equal(model.columns_.sum(axis=0), np.ones(30))
            _check_label_permutations(model.rows_, rows, 3)
            _check_label_permutations(model.columns_, cols, 3)


def test_spectral_biclustering_kluger():
    """Test Kluger methods on a checkerboard dataset."""
    random_state = 0
    for noise in (0, 0.5):
        S, rows, cols = make_checkerboard((30, 30), 3, noise=noise,
                                          random_state=random_state)
        for method in ('scale', 'bistochastic', 'log'):
            for svd_method in ('randomized', 'arpack'):
                model = SpectralBiclustering(n_clusters=(3, 3),
                                             method=method,
                                             svd_method=svd_method,
                                             random_state=random_state)
                model.fit(S)

                assert_equal(model.rows_.shape, (9, 30))
                assert_equal(model.columns_.shape, (9, 30))
                assert_array_equal(model.rows_.sum(axis=0), np.repeat(3, 30))
                assert_array_equal(model.columns_.sum(axis=0), np.repeat(3, 30))
                _check_label_permutations(model.rows_, rows, 3)
                _check_label_permutations(model.columns_, cols, 3)


def _do_scale_test(scaled):
    """Ensures that rows sum to one constant, and columns to another
    constant.

    """
    row_sum = scaled.sum(axis=1)
    col_sum = scaled.sum(axis=0)
    assert_array_almost_equal(row_sum, np.tile(row_sum.mean(), 100),
                              decimal=1)
    assert_array_almost_equal(col_sum, np.tile(col_sum.mean(), 100),
                              decimal=1)


def _do_bistochastic_test(scaled):
    """Ensure that rows and columns sum to the same constant."""
    _do_scale_test(scaled)
    assert_almost_equal(scaled.sum(axis=0).mean(),
                        scaled.sum(axis=1).mean(),
                        decimal=1)


def test_scale_preprocess():
    generator = check_random_state(0)
    x = generator.rand(100, 100)
    scaled, _, _ = _scale_preprocess(x)
    _do_scale_test(scaled)


def test_bistochastic_preprocess():
    generator = check_random_state(0)
    x = generator.rand(100, 100)
    scaled = _bistochastic_preprocess(x)
    _do_bistochastic_test(scaled)


def test_log_preprocess():
    # adding any constant to a log-scaled matrix should make it
    # bistochastic
    generator = check_random_state(0)
    x = generator.rand(100, 100)
    scaled = _log_preprocess(x) + 1
    _do_bistochastic_test(scaled)


def test_fit_best_piecewise():
    vectors = np.array([[0, 0, 0, 1, 1, 1],
                        [2, 2, 2, 3, 3, 3],
                        [0, 1, 2, 3, 4, 5],])
    best = _fit_best_piecewise(vectors, k=2, n_clusters=2,
                               random_state=0, kmeans_kwargs={})
    assert_array_equal(best, vectors[0:2])


def test_project_and_cluster():
    data = np.array([[1, 1, 1,],
                     [1, 1, 1,],
                     [3, 6, 3,],
                     [3, 6, 3,]])
    vectors = np.array([[1, 0,],
                        [0, 1,],
                        [0, 0,]])
    labels = _project_and_cluster(data, vectors, n_clusters=2,
                                  random_state=0, kmeans_kwargs={})
    assert_array_equal(labels, [0, 0, 1, 1])


def test_errors():
    assert_raises(ValueError, SpectralBiclustering, n_clusters=(3, 3, 3))
    assert_raises(ValueError, SpectralBiclustering, method='unknown')
    assert_raises(ValueError, SpectralBiclustering, svd_method='unknown')
    assert_raises(ValueError, SpectralBiclustering, n_components=3, n_best=4)

    model = SpectralBiclustering()
    data = np.arange(27).reshape((3, 3, 3))
    assert_raises(ValueError, model.fit, data)
