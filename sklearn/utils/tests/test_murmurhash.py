# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from sklearn.utils.murmurhash import _murmurhash3_32, murmurhash3_32


def test_mmhash3_int():
    assert _murmurhash3_32(3) == 847579505
    assert _murmurhash3_32(3, seed=0) == 847579505
    assert _murmurhash3_32(3, seed=42) == -1823081949

    assert _murmurhash3_32(3, positive=False) == 847579505
    assert _murmurhash3_32(3, seed=0, positive=False) == 847579505
    assert _murmurhash3_32(3, seed=42, positive=False) == -1823081949

    assert _murmurhash3_32(3, positive=True) == 847579505
    assert _murmurhash3_32(3, seed=0, positive=True) == 847579505
    assert _murmurhash3_32(3, seed=42, positive=True) == 2471885347


def test_mmhash3_int_array():
    rng = np.random.RandomState(42)
    keys = rng.randint(-5342534, 345345, size=3 * 2 * 1).astype(np.int32)
    keys = keys.reshape((3, 2, 1))

    for seed in [0, 42]:
        expected = np.array([_murmurhash3_32(int(k), seed) for k in keys.flat])
        expected = expected.reshape(keys.shape)
        assert_array_equal(_murmurhash3_32(keys, seed), expected)

    for seed in [0, 42]:
        expected = np.array(
            [_murmurhash3_32(k, seed, positive=True) for k in keys.flat]
        )
        expected = expected.reshape(keys.shape)
        assert_array_equal(_murmurhash3_32(keys, seed, positive=True), expected)


def test_mmhash3_bytes():
    assert _murmurhash3_32(b"foo", 0) == -156908512
    assert _murmurhash3_32(b"foo", 42) == -1322301282

    assert _murmurhash3_32(b"foo", 0, positive=True) == 4138058784
    assert _murmurhash3_32(b"foo", 42, positive=True) == 2972666014


def test_mmhash3_unicode():
    assert _murmurhash3_32("foo", 0) == -156908512
    assert _murmurhash3_32("foo", 42) == -1322301282

    assert _murmurhash3_32("foo", 0, positive=True) == 4138058784
    assert _murmurhash3_32("foo", 42, positive=True) == 2972666014


def test_no_collision_on_byte_range():
    previous_hashes = set()
    for i in range(100):
        h = _murmurhash3_32(" " * i, 0)
        assert h not in previous_hashes, "Found collision on growing empty string"


def test_uniform_distribution():
    n_bins, n_samples = 10, 100000
    bins = np.zeros(n_bins, dtype=np.float64)

    for i in range(n_samples):
        bins[_murmurhash3_32(i, positive=True) % n_bins] += 1

    means = bins / n_samples
    expected = np.full(n_bins, 1.0 / n_bins)

    assert_array_almost_equal(means / expected, np.ones(n_bins), 2)


def test_deprecation_warning():
    with pytest.warns(FutureWarning, match="`murmurhash3_32` was deprecated"):
        murmurhash3_32(3)
