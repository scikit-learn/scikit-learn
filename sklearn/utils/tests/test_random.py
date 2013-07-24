from __future__ import division

import numpy as np
from scipy.misc import comb as combinations

from sklearn.utils.random import sample_without_replacement

from sklearn.utils.testing import (
    assert_raises,
    assert_equal,
    assert_true)


###############################################################################
# test custom sampling without replacement algorithm
###############################################################################
def test_invalid_sample_without_replacement_algorithm():
    assert_raises(ValueError, sample_without_replacement, 5, 4, "unknown")


def test_sample_without_replacement_algorithms():
    methods = ("auto", "tracking_selection", "reservoir_sampling", "pool")

    for m in methods:
        def sample_without_replacement_method(n_population, n_samples,
                                              random_state=None):
            return sample_without_replacement(n_population, n_samples,
                                              method=m,
                                              random_state=random_state)

        check_edge_case_of_sample_int(sample_without_replacement_method)
        check_sample_int(sample_without_replacement_method)
        check_sample_int_distribution(sample_without_replacement_method)


def check_edge_case_of_sample_int(sample_without_replacement):

    # n_poluation < n_sample
    assert_raises(ValueError, sample_without_replacement, 0, 1)
    assert_raises(ValueError, sample_without_replacement, 1, 2)

    # n_population == n_samples
    assert_equal(sample_without_replacement(0, 0).shape, (0, ))

    assert_equal(sample_without_replacement(1, 1).shape, (1, ))

    # n_population >= n_samples
    assert_equal(sample_without_replacement(5, 0).shape, (0, ))
    assert_equal(sample_without_replacement(5, 1).shape, (1, ))

    # n_population < 0 or n_samples < 0
    assert_raises(ValueError, sample_without_replacement, -1, 5)
    assert_raises(ValueError, sample_without_replacement, 5, -1)


def check_sample_int(sample_without_replacement):
    # This test is heavily inspired from test_random.py of python-core.
    #
    # For the entire allowable range of 0 <= k <= N, validate that
    # the sample is of the correct length and contains only unique items
    n_population = 100

    for n_samples in range(n_population + 1):
        s = sample_without_replacement(n_population, n_samples)
        assert_equal(len(s), n_samples)
        unique = np.unique(s)
        assert_equal(np.size(unique), n_samples)
        assert_true(np.all(unique < n_population))

    # test edge case n_population == n_samples == 0
    assert_equal(np.size(sample_without_replacement(0, 0)), 0)


def check_sample_int_distribution(sample_without_replacement):
    # This test is heavily inspired from test_random.py of python-core.
    #
    # For the entire allowable range of 0 <= k <= N, validate that
    # sample generates all possible permutations
    n_population = 10

    # a large number of trials prevents false negatives without slowing normal
    # case
    n_trials = 10000

    for n_samples in range(n_population):
        # Counting the number of combinations is not as good as counting the
        # the number of permutations. However, it works with sampling algorithm
        # that does not provide a random permutation of the subset of integer.
        n_expected = combinations(n_population, n_samples, exact=True)

        output = {}
        for i in range(n_trials):
            output[frozenset(sample_without_replacement(n_population,
                                                        n_samples))] = None

            if len(output) == n_expected:
                break
        else:
            raise AssertionError(
                "number of combinations != number of expected (%s != %s)" %
                (len(output), n_expected))
