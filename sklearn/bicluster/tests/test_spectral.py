"""Testing for Spectral Bilustering methods"""

from sklearn.externals.six.moves import cPickle
from sklearn.metrics.pairwise import kernel_metrics

dumps, loads = cPickle.dumps, cPickle.loads

import numpy as np
from scipy import sparse

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater

from sklearn.bicluster.spectral import \
    SpectralBiclustering, \
    scale_preprocess, bistochastic_preprocess, \
    log_preprocess, fit_best_piecewise, \
    project_and_cluster

def test_spectral_biclustering_dhillon():
    S = np.zeros((30, 30))
    S[0:10, 0:10] = 1
    S[10:20, 10:20] = 2
    S[20:30, 20:30] = 3

    model = SpectralBiclustering(random_state=0, n_clusters=3,
                                 method='dhillon')
    model.fit(S)
    assert_equal(model.rows_.shape, (3, 30))
    assert_array_equal(model.rows_.sum(axis=0), np.ones(30))
    assert_array_equal(model.rows_.sum(axis=1), [10] * 3)
    assert_array_equal(model.columns_.sum(axis=0), np.ones(30))
    assert_array_equal(model.columns_.sum(axis=1), [10] * 3)

    model_copy = loads(dumps(model))
    assert_equal(model_copy.n_clusters, model.n_clusters)
    assert_equal(model_copy.method, model.method)


def test_spectral_biclustering_kluger():
    row_vector = [1] * 10 + [3] * 10 + [5] * 10
    col_vector = [2] * 10 + [4] * 10 + [6] * 10
    S = np.outer(row_vector, col_vector)

    for method in ('scale', 'bistochastic', 'log'):
        model = SpectralBiclustering(random_state=0, n_clusters=(3, 3),
                                     method=method)
        model.fit(S)
        assert_equal(model.rows_.shape, (9, 30))
        assert_array_equal(model.rows_.sum(axis=0), 3 * np.ones(30))
        assert_array_equal(model.rows_.sum(axis=1), [10] * 9)

        assert_equal(model.columns_.shape, (9, 30))
        assert_array_equal(model.columns_.sum(axis=0), 3 * np.ones(30))
        assert_array_equal(model.columns_.sum(axis=1), [10] * 9)

        model_copy = loads(dumps(model))
        assert_equal(model_copy.n_clusters, model.n_clusters)
        assert_equal(model_copy.method, model.method)
