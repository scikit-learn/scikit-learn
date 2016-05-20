import numpy as np
from sklearn.decomposition import RobustPCA
from sklearn.utils.testing import assert_array_almost_equal


def _gen_synthetic(n, density, rank):
    r = rank  # Rank
    X = np.random.normal(0, 1/float(n), size=(n, r))
    Y = np.random.normal(0, 1/float(n), size=(n, r))
    L = np.dot(X, Y.T)
    p = density/2
    S = np.random.choice([0, 1, -1], size=(n, n), p=[1 - 2*p, p, p])

    return L, S


def test_rpca():
    rpca = RobustPCA()
    rank = 25
    L, S = _gen_synthetic(500, 0.05, rank)
    X = L + S
    rpca.fit(X)

    assert rpca.n_components_ == rank
    assert_array_almost_equal(L, rpca.low_rank_)
    assert_array_almost_equal(L, rpca.inverse_transform(rpca.transform(X)))
