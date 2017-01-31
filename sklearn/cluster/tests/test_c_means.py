"""Testing for C-means"""

import numpy as np

from sklearn.utils.testing import assert_raises_regex

from sklearn.cluster import CMeans, c_means


def test_c_means_n_init():
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(40, 2))

    # two regression tests on bad n_init argument
    # previous bug: n_init <= 0 threw non-informative TypeError (#3858)
    assert_raises_regex(ValueError, "n_init", CMeans(n_init=0).fit, X)
    assert_raises_regex(ValueError, "n_init", CMeans(n_init=-1).fit, X)