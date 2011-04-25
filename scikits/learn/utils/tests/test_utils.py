import numpy as np

from scikits.learn.utils import make_rng
from nose.tools import assert_raises


def test_make_rng():
    """Check the make_rng utility function behavior"""
    assert make_rng(None) is np.random
    assert make_rng(np.random) is np.random

    rng_42 = np.random.RandomState(42)
    assert make_rng(42).randint(100) == rng_42.randint(100)

    rng_42 = np.random.RandomState(42)
    assert make_rng(rng_42) is rng_42

    rng_42 = np.random.RandomState(42)
    assert make_rng(43).randint(100) != rng_42.randint(100)

    assert_raises(ValueError, make_rng, "some invalid seed")
