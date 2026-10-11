import numpy as np

from sklearn.tree._utils import PytestWeightedFenwickTree


def test_cython_weighted_fenwick_tree(global_random_seed):
    """
    Test Cython's weighted Fenwick tree implementation
    """
    rng = np.random.default_rng(global_random_seed)

    n = 100
    indices = rng.permutation(n)
    y = rng.normal(size=n)
    w = rng.integers(0, 4, size=n)
    y_included_so_far = np.zeros_like(y)
    w_included_so_far = np.zeros_like(w)

    tree = PytestWeightedFenwickTree(n)
    tree.py_reset(n)

    for i in range(n):
        idx = indices[i]
        tree.py_add(idx, y[idx], w[idx])
        y_included_so_far[idx] = y[idx]
        w_included_so_far[idx] = w[idx]

        target = rng.uniform(0, w_included_so_far.sum())
        t_idx_low, t_idx, cw, cwy = tree.py_search(target)

        # check the aggregates are consistent with the returned idx
        assert np.isclose(cw, np.sum(w_included_so_far[:t_idx]))
        assert np.isclose(
            cwy, np.sum(w_included_so_far[:t_idx] * y_included_so_far[:t_idx])
        )

        # check if the cumulative weight is less than or equal to the target
        # depending on t_idx_low and t_idx
        if t_idx_low == t_idx:
            assert cw < target
        else:
            assert cw == target

        # check that if we add the next non-null weight, we are above the target:
        next_weights = w_included_so_far[t_idx:][w_included_so_far[t_idx:] > 0]
        if next_weights.size > 0:
            assert cw + next_weights[0] > target
        # and not below the target for `t_idx_low`:
        next_weights = w_included_so_far[t_idx_low:][w_included_so_far[t_idx_low:] > 0]
        if next_weights.size > 0:
            assert cw + next_weights[0] >= target
