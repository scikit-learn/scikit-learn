# encoding: utf-8
#
# Author: Arnaud Joly
#
# License: BSD Style.
"""
This module contains fast function used to implement random projection
matrix.
"""
cimport cython

import numpy as np
cimport numpy as np
np.import_array()

from sklearn.utils import check_random_state

cpdef _sample_int_check_input(np.int_t n_population, np.int_t n_samples):
    """ Check that input are consistent for sample_int function"""
    if n_population < 0:
        raise ValueError('n_population should be greater than 0, got %s.'
                         % n_population)

    if n_samples > n_population:
        raise ValueError('n_population should be greater or equal than '
                         'n_samples, got n_samples > n_population (%s > %s)'
                         % (n_samples, n_population))

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sample_int_with_set_based_method(np.int_t n_population,
                                       np.int_t n_samples,
                                       random_state=None):
    """Sample integers without replacement.

    Select n_samples integers from the set [0, n_population) without
    replacement.

    Parameters
    ----------
    n_population: int,
        The size of the set to sample from.

    n_samples: int,
        The number of integer to sample.

    random_state: int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    result : array of size (n_samples, )
        The sampled subsets of integer.
    """
    _sample_int_check_input(n_population, n_samples)

    cdef np.int_t i
    cdef np.int_t j
    cdef np.ndarray[np.int_t, ndim=1] result = np.empty((n_samples, ),
                                                        dtype=np.int)

    rng = check_random_state(random_state)
    rng_randint = rng.randint

    # The following line of code are heavily inspired from python core,
    # more precisely of random.sample.
    cdef set selected = set()
    selected_add = selected.add

    for i in range(n_samples):
        j = rng_randint(n_population)
        while j in selected:
            j = rng_randint(n_population)
        selected_add(j)
        result[i] = j

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sample_int_with_reservoir_sampling(np.int_t n_population,
                                         np.int_t n_samples,
                                         random_state=None):
    """Sample integers without replacement.

    Select n_samples integers from the set [0, n_population) without
    replacement.

    Time complexity of O(n_population)
    Space complexity of O(n_samples)


    Parameters
    ----------
    n_population: int,
        The size of the set to sample from.

    n_samples: int,
         The number of integer to sample.

    random_state: int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    result : array of size (n_samples, )
        The sampled subsets of integer. The order of the items is not
        necessarily random. Use a random permutation of the array if the order
        of the items to be randomized.
    """
    _sample_int_check_input(n_population, n_samples)

    cdef np.int_t i
    cdef np.int_t j
    cdef np.ndarray[np.int_t, ndim=1] result = np.empty((n_samples, ),
                                                        dtype=np.int)

    rng = check_random_state(random_state)
    rng_randint = rng.randint

    # This cython implementation  is based on the one of Robert Kern:
    # http://mail.scipy.org/pipermail/numpy-discussion/2010-December/
    # 054289.html
    #
    for i in range(n_samples):
        result[i] = i

    for i from n_samples <= i < n_population:
        j = rng_randint(0, i + 1)
        if j < n_samples:
            result[j] = i

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sample_int(np.int_t n_population, np.int_t n_samples, method="auto",
                 random_state=None):
    """Sample integers without replacement.

    Select n_samples integers from the set [0, n_population) without
    replacement.


    Parameters
    ----------
    n_population: int,
        The size of the set to sample from.

    n_samples: int,
        The number of integer to sample.

    random_state: int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    method: "auto", "set_based" or "reservoir_sampling"
        If method == "auto", an algorithm is automatically selected.
        If method =="set_based", a set based implementation is used which is
        suitable for `n_samples` <<< `n_population`.
        If method == "reservoir_sampling", a reservoir samplingn algorithm
        which is suitable for high memory constraint or when
        `n_samples` ~ `n_population`

    Returns
    -------
    result : array of size (n_samples, )
        The sampled subsets of integer.
    """
    _sample_int_check_input(n_population, n_samples)

    all_methods = ("auto", "set_based", "reservoir_sampling")

    if method == "auto" or method == "set_based":
        return sample_int_with_set_based_method(n_population, n_samples,
                                                random_state)

    elif method == "reservoir_sampling":
        return sample_int_with_reservoir_sampling(n_population, n_samples,
                                                  random_state)
    else:
        raise ValueError('Expected a method name in %s, got %s. '
                         % (all_methods, method))


