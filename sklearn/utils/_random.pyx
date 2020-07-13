# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Arnaud Joly
#
# License: BSD 3 clause
"""
Random utility function
=======================
This module complements missing features of ``numpy.random``.

The module contains:
    * Several algorithms to sample integers without replacement.
    * Fast rand_r alternative based on xor shifts
"""
cimport cython

import numpy as np
cimport numpy as np
np.import_array()

from . import check_random_state

cdef UINT32_t DEFAULT_SEED = 1


cpdef _sample_without_replacement_check_input(np.int_t n_population,
                                              np.int_t n_samples):
    """ Check that input are consistent for sample_without_replacement"""
    if n_population < 0:
        raise ValueError('n_population should be greater than 0, got %s.'
                         % n_population)

    if n_samples > n_population:
        raise ValueError('n_population should be greater or equal than '
                         'n_samples, got n_samples > n_population (%s > %s)'
                         % (n_samples, n_population))


cpdef _sample_without_replacement_with_tracking_selection(
        np.int_t n_population,
        np.int_t n_samples,
        random_state=None):
    r"""Sample integers without replacement.

    Select n_samples integers from the set [0, n_population) without
    replacement.

    Time complexity:
        - Worst-case: unbounded
        - Average-case:
            O(O(np.random.randint) * \sum_{i=1}^n_samples 1 /
                                              (1 - i / n_population)))
            <= O(O(np.random.randint) *
                   n_population * ln((n_population - 2)
                                     /(n_population - 1 - n_samples)))
            <= O(O(np.random.randint) *
                 n_population * 1 / (1 - n_samples / n_population))

    Space complexity of O(n_samples) in a python set.


    Parameters
    ----------
    n_population : int,
        The size of the set to sample from.

    n_samples : int,
        The number of integer to sample.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    out : array of size (n_samples, )
        The sampled subsets of integer.
    """
    _sample_without_replacement_check_input(n_population, n_samples)

    cdef np.int_t i
    cdef np.int_t j
    cdef np.ndarray[np.int_t, ndim=1] out = np.empty((n_samples, ), dtype=int)

    rng = check_random_state(random_state)
    rng_randint = rng.randint

    # The following line of code are heavily inspired from python core,
    # more precisely of random.sample.
    cdef set selected = set()

    for i in range(n_samples):
        j = rng_randint(n_population)
        while j in selected:
            j = rng_randint(n_population)
        selected.add(j)
        out[i] = j

    return out


cpdef _sample_without_replacement_with_pool(np.int_t n_population,
                                            np.int_t n_samples,
                                            random_state=None):
    """Sample integers without replacement.

    Select n_samples integers from the set [0, n_population) without
    replacement.

    Time complexity: O(n_population +  O(np.random.randint) * n_samples)

    Space complexity of O(n_population + n_samples).


    Parameters
    ----------
    n_population : int,
        The size of the set to sample from.

    n_samples : int,
        The number of integer to sample.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    out : array of size (n_samples, )
        The sampled subsets of integer.
    """
    _sample_without_replacement_check_input(n_population, n_samples)

    cdef np.int_t i
    cdef np.int_t j
    cdef np.ndarray[np.int_t, ndim=1] out = np.empty((n_samples, ), dtype=int)

    cdef np.ndarray[np.int_t, ndim=1] pool = np.empty((n_population, ),
                                                      dtype=int)

    rng = check_random_state(random_state)
    rng_randint = rng.randint

    # Initialize the pool
    for i in range(n_population):
        pool[i] = i

    # The following line of code are heavily inspired from python core,
    # more precisely of random.sample.
    for i in range(n_samples):
        j = rng_randint(n_population - i)  # invariant: non-selected at [0,n-i)
        out[i] = pool[j]
        pool[j] = pool[n_population - i - 1]  # move non-selected item into
                                              # vacancy

    return out


cpdef _sample_without_replacement_with_reservoir_sampling(
    np.int_t n_population,
    np.int_t n_samples,
    random_state=None):
    """Sample integers without replacement.

    Select n_samples integers from the set [0, n_population) without
    replacement.

    Time complexity of
        O((n_population - n_samples) * O(np.random.randint) + n_samples)
    Space complexity of O(n_samples)


    Parameters
    ----------
    n_population : int,
        The size of the set to sample from.

    n_samples : int,
         The number of integer to sample.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    out : array of size (n_samples, )
        The sampled subsets of integer. The order of the items is not
        necessarily random. Use a random permutation of the array if the order
        of the items has to be randomized.
    """
    _sample_without_replacement_check_input(n_population, n_samples)

    cdef np.int_t i
    cdef np.int_t j
    cdef np.ndarray[np.int_t, ndim=1] out = np.empty((n_samples, ), dtype=int)

    rng = check_random_state(random_state)
    rng_randint = rng.randint

    # This cython implementation is based on the one of Robert Kern:
    # http://mail.scipy.org/pipermail/numpy-discussion/2010-December/
    # 054289.html
    #
    for i in range(n_samples):
        out[i] = i

    for i from n_samples <= i < n_population:
        j = rng_randint(0, i + 1)
        if j < n_samples:
            out[j] = i

    return out


cpdef sample_without_replacement(np.int_t n_population,
                                 np.int_t n_samples,
                                 method="auto",
                                 random_state=None):
    """Sample integers without replacement.

    Select n_samples integers from the set [0, n_population) without
    replacement.


    Parameters
    ----------
    n_population : int,
        The size of the set to sample from.

    n_samples : int,
        The number of integer to sample.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    method : "auto", "tracking_selection", "reservoir_sampling" or "pool"
        If method == "auto", the ratio of n_samples / n_population is used
        to determine which algorithm to use:
        If ratio is between 0 and 0.01, tracking selection is used.
        If ratio is between 0.01 and 0.99, numpy.random.permutation is used.
        If ratio is greater than 0.99, reservoir sampling is used.
        The order of the selected integers is undefined. If a random order is
        desired, the selected subset should be shuffled.

        If method =="tracking_selection", a set based implementation is used
        which is suitable for `n_samples` <<< `n_population`.

        If method == "reservoir_sampling", a reservoir sampling algorithm is
        used which is suitable for high memory constraint or when
        O(`n_samples`) ~ O(`n_population`).
        The order of the selected integers is undefined. If a random order is
        desired, the selected subset should be shuffled.

        If method == "pool", a pool based algorithm is particularly fast, even
        faster than the tracking selection method. Hovewer, a vector containing
        the entire population has to be initialized.
        If n_samples ~ n_population, the reservoir sampling method is faster.

    Returns
    -------
    out : array of size (n_samples, )
        The sampled subsets of integer. The subset of selected integer might
        not be randomized, see the method argument.
    """
    _sample_without_replacement_check_input(n_population, n_samples)

    all_methods = ("auto", "tracking_selection", "reservoir_sampling", "pool")

    ratio = n_samples / n_population if n_population != 0.0 else 1.0

    # Check ratio and use permutation unless ratio < 0.01 or ratio > 0.99
    if method == "auto" and ratio > 0.01 and ratio < 0.99:
        rng = check_random_state(random_state)
        return rng.permutation(n_population)[:n_samples]

    if method == "auto" or method == "tracking_selection":
        # TODO the pool based method can also be used.
        #      however, it requires special benchmark to take into account
        #      the memory requirement of the array vs the set.

        # The value 0.2 has been determined through benchmarking.
        if ratio < 0.2:
            return _sample_without_replacement_with_tracking_selection(
                n_population, n_samples, random_state)
        else:
            return _sample_without_replacement_with_reservoir_sampling(
                n_population, n_samples, random_state)

    elif method == "reservoir_sampling":
        return _sample_without_replacement_with_reservoir_sampling(
            n_population, n_samples, random_state)

    elif method == "pool":
        return _sample_without_replacement_with_pool(n_population, n_samples,
                                                     random_state)
    else:
        raise ValueError('Expected a method name in %s, got %s. '
                         % (all_methods, method))


def _our_rand_r_py(seed):
    """Python utils to test the our_rand_r function"""
    cdef UINT32_t my_seed = seed
    return our_rand_r(&my_seed)
