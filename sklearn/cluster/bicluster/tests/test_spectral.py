"""Testing for Spectral Biclustering methods"""

import numpy as np
from scipy.sparse import csr_matrix, issparse

from sklearn.grid_search import ParameterGrid

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises

from sklearn.cluster.bicluster import SpectralCoclustering
from sklearn.cluster.bicluster import SpectralBiclustering
from sklearn.cluster.bicluster.spectral import _scale_normalize
from sklearn.cluster.bicluster.spectral import _bistochastic_normalize
from sklearn.cluster.bicluster.spectral import _log_normalize

from sklearn.metrics.cluster.bicluster import consensus_score

from sklearn.datasets import make_biclusters, make_checkerboard


def test_spectral_coclustering():
    """Test Dhillon's Spectral CoClustering on a simple problem."""
    param_grid = {'svd_method': ['randomized', 'arpack'],
                  'n_svd_vecs': [None, 20],
                  'mini_batch': [False, True],
                  'init': ['k-means++'],
                  'n_init': [10],
                  'n_jobs': [1]}
    random_state = 0
    S, rows, cols = make_biclusters((30, 30), 3, noise=0.5,
                                    random_state=random_state)
    S -= S.min()  # needs to be nonnegative before making it sparse
    S = np.where(S < 1, 0, S)  # threshold some values
    for mat in (S, csr_matrix(S)):
        for kwargs in ParameterGrid(param_grid):
            model = SpectralCoclustering(n_clusters=3,
                                         random_state=random_state,
                                         **kwargs)
            model.fit(mat)

            assert_equal(model.rows_.shape, (3, 30))
            assert_array_equal(model.rows_.sum(axis=0), np.ones(30))
            assert_array_equal(model.columns_.sum(axis=0), np.ones(30))
            assert_equal(consensus_score(model.biclusters_,
                                         (rows, cols)), 1)


def test_spectral_biclustering():
    """Test Kluger methods on a checkerboard dataset."""
    param_grid = {'method': ['scale', 'bistochastic', 'log'],
                  'svd_method': ['randomized', 'arpack'],
                  'n_svd_vecs': [None, 20],
                  'mini_batch': [False, True],
                  'init': ['k-means++'],
                  'n_init': [10],
                  'n_jobs': [1]}
    random_state = 0
    S, rows, cols = make_checkerboard((30, 30), 3, noise=0.5,
                                      random_state=random_state)
    for mat in (S, csr_matrix(S)):
        for kwargs in ParameterGrid(param_grid):
            model = SpectralBiclustering(n_clusters=3,
                                         random_state=random_state,
                                         **kwargs)

            if issparse(mat) and kwargs['method'] == 'log':
                # cannot take log of sparse matrix
                assert_raises(ValueError, model.fit, mat)
                continue
            else:
                model.fit(mat)

            assert_equal(model.rows_.shape, (9, 30))
            assert_equal(model.columns_.shape, (9, 30))
            assert_array_equal(model.rows_.sum(axis=0),
                               np.repeat(3, 30))
            assert_array_equal(model.columns_.sum(axis=0),
                               np.repeat(3, 30))
            assert_equal(consensus_score(model.biclusters_,
                                         (rows, cols)), 1)


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


def test_scale_normalize():
    generator = np.random.RandomState(0)
    X = generator.rand(100, 100)
    for mat in (X, csr_matrix(X)):
        scaled, _, _ = _scale_normalize(mat)
        _do_scale_test(scaled)
        if issparse(mat):
            assert issparse(scaled)


def test_bistochastic_normalize():
    generator = np.random.RandomState(0)
    X = generator.rand(100, 100)
    for mat in (X, csr_matrix(X)):
        scaled = _bistochastic_normalize(mat)
        _do_bistochastic_test(scaled)
        if issparse(mat):
            assert issparse(scaled)


def test_log_normalize():
    # adding any constant to a log-scaled matrix should make it
    # bistochastic
    generator = np.random.RandomState(0)
    mat = generator.rand(100, 100)
    scaled = _log_normalize(mat) + 1
    _do_bistochastic_test(scaled)


def test_fit_best_piecewise():
    model = SpectralBiclustering(random_state=0)
    vectors = np.array([[0, 0, 0, 1, 1, 1],
                        [2, 2, 2, 3, 3, 3],
                        [0, 1, 2, 3, 4, 5]])
    best = model._fit_best_piecewise(vectors, n_best=2, n_clusters=2)
    assert_array_equal(best, vectors[:2])


def test_project_and_cluster():
    model = SpectralBiclustering(random_state=0)
    data = np.array([[1, 1, 1],
                     [1, 1, 1],
                     [3, 6, 3],
                     [3, 6, 3]])
    vectors = np.array([[1, 0],
                        [0, 1],
                        [0, 0]])
    for mat in (data, csr_matrix(data)):
        labels = model._project_and_cluster(data, vectors,
                                            n_clusters=2)
        assert_array_equal(labels, [0, 0, 1, 1])


def test_perfect_checkerboard():
    model = SpectralBiclustering(3, svd_method="arpack", random_state=0)

    S, rows, cols = make_checkerboard((30, 30), 3, noise=0,
                                      random_state=0)
    model.fit(S)
    assert_equal(consensus_score(model.biclusters_,
                                 (rows, cols)), 1)

    S, rows, cols = make_checkerboard((40, 30), 3, noise=0,
                                      random_state=0)
    model.fit(S)
    assert_equal(consensus_score(model.biclusters_,
                                 (rows, cols)), 1)

    S, rows, cols = make_checkerboard((30, 40), 3, noise=0,
                                      random_state=0)
    model.fit(S)
    assert_equal(consensus_score(model.biclusters_,
                                 (rows, cols)), 1)


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

    model = SpectralBiclustering(n_components=0)
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering(n_best=0)
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering(n_components=3, n_best=4)
    assert_raises(ValueError, model.fit, data)

    model = SpectralBiclustering()
    data = np.arange(27).reshape((3, 3, 3))
    assert_raises(ValueError, model.fit, data)
