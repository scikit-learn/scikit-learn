# Author: Hamzeh Alsalhi <ha258@cornell.edu>
#
# License: BSD 3 clause
from __future__ import division
import numpy as np
import scipy.sparse as sp
import operator
import array

from sklearn.utils import check_random_state

from ._random import sample_without_replacement

__all__ = ['sample_without_replacement', 'choice']


# This is a backport of np.random.choice from numpy 1.7
# The function can be removed when we bump the requirements to >=1.7
def choice(a, size=None, replace=True, p=None, random_state=None):
    """
    choice(a, size=None, replace=True, p=None)

    Generates a random sample from a given 1-D array

    .. versionadded:: 1.7.0

    Parameters
    -----------
    a : 1-D array-like or int
        If an ndarray, a random sample is generated from its elements.
        If an int, the random sample is generated as if a was np.arange(n)

    size : int or tuple of ints, optional
        Output shape. Default is None, in which case a single value is
        returned.

    replace : boolean, optional
        Whether the sample is with or without replacement.

    p : 1-D array-like, optional
        The probabilities associated with each entry in a.
        If not given the sample assumes a uniform distribtion over all
        entries in a.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.


    Returns
    --------
    samples : 1-D ndarray, shape (size,)
    The generated random samples

    Raises
    -------
    ValueError
    If a is an int and less than zero, if a or p are not 1-dimensional,
    if a is an array-like of size 0, if p is not a vector of
    probabilities, if a and p have different lengths, or if
    replace=False and the sample size is greater than the population
    size

    See Also
    ---------
    randint, shuffle, permutation

    Examples
    ---------
    Generate a uniform random sample from np.arange(5) of size 3:

    >>> np.random.choice(5, 3)  # doctest: +SKIP
    array([0, 3, 4])
    >>> #This is equivalent to np.random.randint(0,5,3)

    Generate a non-uniform random sample from np.arange(5) of size 3:

    >>> np.random.choice(5, 3, p=[0.1, 0, 0.3, 0.6, 0])  # doctest: +SKIP
    array([3, 3, 0])

    Generate a uniform random sample from np.arange(5) of size 3 without
    replacement:

    >>> np.random.choice(5, 3, replace=False)  # doctest: +SKIP
    array([3,1,0])
    >>> #This is equivalent to np.random.shuffle(np.arange(5))[:3]

    Generate a non-uniform random sample from np.arange(5) of size
    3 without replacement:

    >>> np.random.choice(5, 3, replace=False, p=[0.1, 0, 0.3, 0.6, 0])
    ... # doctest: +SKIP
    array([2, 3, 0])

    Any of the above can be repeated with an arbitrary array-like
    instead of just integers. For instance:

    >>> aa_milne_arr = ['pooh', 'rabbit', 'piglet', 'Christopher']
    >>> np.random.choice(aa_milne_arr, 5, p=[0.5, 0.1, 0.1, 0.3])
    ... # doctest: +SKIP
    array(['pooh', 'pooh', 'pooh', 'Christopher', 'piglet'],
    dtype='|S11')

    """
    random_state = check_random_state(random_state)

    # Format and Verify input
    a = np.array(a, copy=False)
    if a.ndim == 0:
        try:
            # __index__ must return an integer by python rules.
            pop_size = operator.index(a.item())
        except TypeError:
            raise ValueError("a must be 1-dimensional or an integer")
        if pop_size <= 0:
            raise ValueError("a must be greater than 0")
    elif a.ndim != 1:
        raise ValueError("a must be 1-dimensional")
    else:
        pop_size = a.shape[0]
        if pop_size is 0:
            raise ValueError("a must be non-empty")

    if None != p:
        p = np.array(p, dtype=np.double, ndmin=1, copy=False)
        if p.ndim != 1:
            raise ValueError("p must be 1-dimensional")
        if p.size != pop_size:
            raise ValueError("a and p must have same size")
        if np.any(p < 0):
            raise ValueError("probabilities are not non-negative")
        if not np.allclose(p.sum(), 1):
            raise ValueError("probabilities do not sum to 1")

    shape = size
    if shape is not None:
        size = np.prod(shape, dtype=np.intp)
    else:
        size = 1

    # Actual sampling
    if replace:
        if None != p:
            cdf = p.cumsum()
            cdf /= cdf[-1]
            uniform_samples = random_state.random_sample(shape)
            idx = cdf.searchsorted(uniform_samples, side='right')
            # searchsorted returns a scalar
            idx = np.array(idx, copy=False)
        else:
            idx = random_state.randint(0, pop_size, size=shape)
    else:
        if size > pop_size:
            raise ValueError("Cannot take a larger sample than "
                             "population when 'replace=False'")

        if None != p:
            if np.sum(p > 0) < size:
                raise ValueError("Fewer non-zero entries in p than size")
            n_uniq = 0
            p = p.copy()
            found = np.zeros(shape, dtype=np.int)
            flat_found = found.ravel()
            while n_uniq < size:
                x = random_state.rand(size - n_uniq)
                if n_uniq > 0:
                    p[flat_found[0:n_uniq]] = 0
                cdf = np.cumsum(p)
                cdf /= cdf[-1]
                new = cdf.searchsorted(x, side='right')
                _, unique_indices = np.unique(new, return_index=True)
                unique_indices.sort()
                new = new.take(unique_indices)
                flat_found[n_uniq:n_uniq + new.size] = new
                n_uniq += new.size
            idx = found
        else:
            idx = random_state.permutation(pop_size)[:size]
            if shape is not None:
                idx.shape = shape

    if shape is None and isinstance(idx, np.ndarray):
        # In most cases a scalar will have been made an array
        idx = idx.item(0)

    # Use samples as indices for a if a is array-like
    if a.ndim == 0:
        return idx

    if shape is not None and idx.ndim == 0:
        # If size == () then the user requested a 0-d array as opposed to
        # a scalar object when size is None. However a[idx] is always a
        # scalar and not an array. So this makes sure the result is an
        # array, taking into account that np.array(item) may not work
        # for object arrays.
        res = np.empty((), dtype=a.dtype)
        res[()] = a[idx]
        return res

    return a[idx]


def random_choice_csc(n_samples, classes, class_probability=None,
                      random_state=None):
    """Generate a sparse random matrix given column class distributions

    Parameters
    ----------
    n_samples : int,
        Number of samples to draw in each column.

    classes : list of size n_outputs of arrays of size (n_classes,)
        List of classes for each column.

    class_probability : list of size n_outputs of arrays of size (n_classes,)
        Optional (default=None). Class distribution of each column. If None the
        uniform distribution is assumed.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    random_matrix : sparse csc matrix of size (n_samples, n_outputs)

    """
    data = array.array('i')
    indices = array.array('i')
    indptr = array.array('i', [0])

    for j in range(len(classes)):
        classes[j] = np.asarray(classes[j])
        if classes[j].dtype.kind != 'i':
            raise ValueError("class dtype %s is not supported" %
                             classes[j].dtype)
        classes[j] = classes[j].astype(int)

        # use uniform distribution if no class_probability is given
        if class_probability is None:
            class_prob_j = np.empty(shape=classes[j].shape[0])
            class_prob_j.fill(1 / classes[j].shape[0])
        else:
            class_prob_j = np.asarray(class_probability[j])

        if np.sum(class_prob_j) != 1.0:
            raise ValueError("Probability array at index {0} does not sum to "
                             "one".format(j))

        if class_prob_j.shape[0] != classes[j].shape[0]:
            raise ValueError("classes[{0}] (length {1}) and "
                             "class_probability[{0}] (length {2}) have "
                             "different length.".format(j,
                                                        classes[j].shape[0],
                                                        class_prob_j.shape[0]))

        # If 0 is not present in the classes insert it with a probability 0.0
        if 0 not in classes[j]:
            classes[j] = np.insert(classes[j], 0, 0)
            class_prob_j = np.insert(class_prob_j, 0, 0.0)

        # If there are nonzero classes choose randomly using class_probability
        if classes[j].shape[0] > 1:
            p_nonzero = 1 - class_prob_j[classes[j] == 0]
            nnz = int(n_samples * p_nonzero)
            ind_sample = sample_without_replacement(n_population=n_samples,
                                                    n_samples=nnz,
                                                    random_state=random_state)
            indices.extend(ind_sample)

            # Normalize probabilites for the nonzero elements
            classes_j_nonzero = classes[j] != 0
            class_probability_nz = class_prob_j[classes_j_nonzero]
            class_probability_nz_norm = (class_probability_nz /
                                         np.sum(class_probability_nz))
            classes_ind = np.searchsorted(class_probability_nz_norm.cumsum(),
                                          np.random.rand(nnz))
            data.extend(classes[j][classes_j_nonzero][classes_ind])
        indptr.append(len(indices))

    return sp.csc_matrix((data, indices, indptr),
                         (n_samples, len(classes)),
                         dtype=int)
