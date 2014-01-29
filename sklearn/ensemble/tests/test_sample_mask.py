"""
Testing for the gradient boosting module (sklearn.ensemble.gradient_boosting).
"""

import numpy as np

from sklearn.ensemble._gradient_boosting import _ranked_random_sample_mask

rng = np.random.RandomState(0)


def test_ranked_sample_mask():
    """Check sampling on a groups."""
    n_groups = 10
    for t in [np.int32, np.int64]:
        group = np.round(np.arange(0, n_groups - 1, 0.2)).astype(t)
        n_total_samples = len(group)
        n_uniq_group = len(np.unique(group))
        n_inbag = 5

        mask = _ranked_random_sample_mask(n_total_samples, n_inbag,
                                          group, n_uniq_group, rng)

        inbag = np.unique(group[mask])
        oob = np.unique(group[~mask])
        assert(len(inbag) == n_inbag)
        assert(set(inbag).intersection(set(oob)) == set())


if __name__ == "__main__":
    import nose
    nose.runmodule()
