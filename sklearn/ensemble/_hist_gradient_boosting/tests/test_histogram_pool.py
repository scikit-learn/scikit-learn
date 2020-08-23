from itertools import product

from numpy.testing import assert_allclose

from sklearn.ensemble._hist_gradient_boosting._histogram_pool import (
    HistogramsPool
)


def test_histograms_pool():
    # simple check how HistogramsPool manages state
    n_features, n_bins = 20, 5
    pool = HistogramsPool(n_features=n_features, n_bins=n_bins)

    histograms1_ref = pool.get()
    histograms1 = histograms1_ref()
    assert histograms1.shape == (n_features, n_bins)
    histograms1[:] = 1

    assert pool.used_pool == [histograms1]
    assert pool.avaliable_pool == []

    histograms2_ref = pool.get()
    histograms2 = histograms2_ref()
    histograms2[:] = 2

    assert histograms2.shape == (n_features, n_bins)
    assert pool.used_pool == [histograms1, histograms2]
    assert pool.avaliable_pool == []

    # when pool is reset histograms in the used pool is moved to the
    # avaliable pool
    pool.reset()
    assert pool.avaliable_pool == [histograms1, histograms2]
    assert pool.used_pool == []

    # histograms are reset to zero
    keys = ['sum_gradients', 'sum_hessians', 'count']
    for histograms, key in product(pool.avaliable_pool, keys):
        assert_allclose(histograms[key], 0)

    histograms3_ref = pool.get()
    histograms3 = histograms3_ref()
    assert histograms3.shape == (n_features, n_bins)

    # only histograms1 is in the avaliable pool
    assert pool.avaliable_pool == [histograms1]
    assert pool.used_pool == [histograms3]
