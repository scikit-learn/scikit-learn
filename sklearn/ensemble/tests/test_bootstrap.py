"""
Testing for the utility function _get_n_samples_bootstrap
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

import numpy as np
import pytest

from sklearn.ensemble._bootstrap import _get_n_samples_bootstrap


def test_get_n_samples_bootstrap():
    # max_samples=None returns n_samples
    n_samples, max_samples, sample_weight = 10, None, "not_used"
    assert _get_n_samples_bootstrap(n_samples, max_samples, sample_weight) == n_samples

    # max_samples:int returns max_samples
    n_samples, max_samples, sample_weight = 10, 5, "not_used"
    assert (
        _get_n_samples_bootstrap(n_samples, max_samples, sample_weight) == max_samples
    )

    # cases where n_samples_bootstrap is small and should raise a warning
    warning_msg = ".+the number of samples.+low number.+max_samples.+as an integer"
    n_samples, max_samples, sample_weight = 10, 0.66, None
    with pytest.warns(UserWarning, match=warning_msg):
        assert _get_n_samples_bootstrap(n_samples, max_samples, sample_weight) == int(
            max_samples * n_samples
        )

    n_samples, max_samples, sample_weight = 10, 0.01, None
    with pytest.warns(UserWarning, match=warning_msg):
        assert _get_n_samples_bootstrap(n_samples, max_samples, sample_weight) == 1

    warning_msg_with_weights = (
        ".+the total sum of sample weights.+low number.+max_samples.+as an integer"
    )
    rng = np.random.default_rng(0)
    n_samples, max_samples, sample_weight = 10, 0.8, rng.uniform(size=10)
    with pytest.warns(UserWarning, match=warning_msg_with_weights):
        assert _get_n_samples_bootstrap(n_samples, max_samples, sample_weight) == int(
            max_samples * sample_weight.sum()
        )

    # cases where n_samples_bootstrap is big enough and shouldn't raise a warning
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        n_samples, max_samples, sample_weight = 100, 30, None
        assert (
            _get_n_samples_bootstrap(n_samples, max_samples, sample_weight)
            == max_samples
        )
        n_samples, max_samples, sample_weight = 100, 0.5, rng.uniform(size=100)
        assert _get_n_samples_bootstrap(n_samples, max_samples, sample_weight) == int(
            max_samples * sample_weight.sum()
        )


@pytest.mark.parametrize("max_samples", [None, 1, 5, 1000, 0.1, 1.0, 1.5])
def test_n_samples_bootstrap_repeated_weighted_equivalence(max_samples):
    # weighted dataset
    n_samples = 100
    rng = np.random.RandomState(0)
    sample_weight = rng.randint(2, 5, n_samples)
    # repeated dataset
    n_samples_repeated = sample_weight.sum()

    n_bootstrap_weighted = _get_n_samples_bootstrap(
        n_samples, max_samples, sample_weight
    )
    n_bootstrap_repeated = _get_n_samples_bootstrap(
        n_samples_repeated, max_samples, None
    )
    if max_samples is None:
        assert n_bootstrap_weighted != n_bootstrap_repeated
    else:
        assert n_bootstrap_weighted == n_bootstrap_repeated
