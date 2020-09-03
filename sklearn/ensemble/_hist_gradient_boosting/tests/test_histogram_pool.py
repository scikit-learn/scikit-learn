from sklearn.ensemble._hist_gradient_boosting._histogram_pool import (
    HistogramPool
)
import pytest


def test_histograms_pool():
    # simple check how HistogramPool manages state
    n_features, n_bins = 20, 5
    pool = HistogramPool(n_features=n_features, n_bins=n_bins)

    histograms1 = pool.get()
    assert histograms1.shape == (n_features, n_bins)

    assert pool.used_pool == [histograms1]
    assert pool.available_pool == []

    histograms2 = pool.get()
    assert histograms2.shape == (n_features, n_bins)
    assert pool.used_pool == [histograms1, histograms2]
    assert pool.available_pool == []

    pool.release(histograms1)
    assert pool.used_pool == [histograms2]
    assert pool.available_pool == [histograms1]

    # Cannot release an already released histogram
    with pytest.raises(ValueError):
        pool.release(histograms1)

    # when pool is reset histograms in the used pool is moved to the
    # avaliable pool
    pool.reset()
    assert pool.available_pool == [histograms1, histograms2]
    assert pool.used_pool == []

    histograms3 = pool.get()
    assert histograms3.shape == (n_features, n_bins)

    # only histograms1 is in the avaliable pool
    assert pool.available_pool == [histograms1]
    assert histograms3 is histograms2
    assert pool.used_pool == [histograms3]
