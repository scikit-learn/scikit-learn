import numpy as np

from sklearn.tree._utils import PytestWeightedFenwickTree


# @pytest.mark.parametrize("min_heap", [True, False])
def test_cython_weighted_fenwick_tree():
    """
    Test Cython's weighted Fenwick tree implementation
    """
    rng = np.random.default_rng()

    n = 100
    indices = rng.permutation(n)
    y = rng.normal(size=n)
    w = rng.integers(1, 4, size=n)
    y_sorted = np.zeros_like(y)
    w_sorted = np.zeros_like(w)

    tree = PytestWeightedFenwickTree(n)
    tree.py_reset(n)

    for idx in indices:
        tree.py_add(idx, y[idx], w[idx])
        y_sorted[idx] = y[idx]
        w_sorted[idx] = w[idx]
        t = rng.uniform(0, w_sorted.sum())
        t_idx, cw, cwy = tree.py_search(t)
        assert np.isclose(cw, w_sorted[:t_idx].sum())
        assert np.isclose(cwy, (w_sorted[:t_idx] * y_sorted[:t_idx]).sum())
        assert cw <= t
