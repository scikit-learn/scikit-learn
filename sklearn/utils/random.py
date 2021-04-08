# Author: Hamzeh Alsalhi <ha258@cornell.edu>
#
# License: BSD 3 clause
import numpy as np
import scipy.sparse as sp
import array

from . import check_random_state
from ._random import sample_without_replacement

__all__ = ['sample_without_replacement']


def _random_choice_csc(n_samples, classes, class_probability=None,
                       random_state=None):
    """Generate a sparse random matrix given column class distributions

    Parameters
    ----------
    n_samples : int,
        Number of samples to draw in each column.

    classes : list of size n_outputs of arrays of size (n_classes,)
        List of classes for each column.

    class_probability : list of size n_outputs of arrays of \
        shape (n_classes,), default=None
        Class distribution of each column. If None, uniform distribution is
        assumed.

    random_state : int, RandomState instance or None, default=None
        Controls the randomness of the sampled classes.
        See :term:`Glossary <random_state>`.

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
        classes[j] = classes[j].astype(np.int64, copy=False)

        # use uniform distribution if no class_probability is given
        if class_probability is None:
            class_prob_j = np.empty(shape=classes[j].shape[0])
            class_prob_j.fill(1 / classes[j].shape[0])
        else:
            class_prob_j = np.asarray(class_probability[j])

        if not np.isclose(np.sum(class_prob_j), 1.0):
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
        rng = check_random_state(random_state)
        if classes[j].shape[0] > 1:
            p_nonzero = 1 - class_prob_j[classes[j] == 0]
            nnz = int(n_samples * p_nonzero)
            ind_sample = sample_without_replacement(n_population=n_samples,
                                                    n_samples=nnz,
                                                    random_state=random_state)
            indices.extend(ind_sample)

            # Normalize probabilities for the nonzero elements
            classes_j_nonzero = classes[j] != 0
            class_probability_nz = class_prob_j[classes_j_nonzero]
            class_probability_nz_norm = (class_probability_nz /
                                         np.sum(class_probability_nz))
            classes_ind = np.searchsorted(class_probability_nz_norm.cumsum(),
                                          rng.rand(nnz))
            data.extend(classes[j][classes_j_nonzero][classes_ind])
        indptr.append(len(indices))

    return sp.csc_matrix((data, indices, indptr),
                         (n_samples, len(classes)),
                         dtype=int)


def sample_erdos_renyi_gnm(n, m, samples=1, random_state=None, method="auto",
                           return_as="edge_array", dtype="int64"):
    """Sample the `G(n, m)`-model of Erdos–Renyi random graphs.

    The problem is solved by the denoting the edges by `0, ..., n(n-1)//2 - 1`,
    sampling those and mapping the result onto the respective indices of the
    adjacency matrix.


    Parameters
    ----------
    n : int
        Number of vertices of the graph.

    m : int
        Number of distinct edges to be in each sample.
        Must fulfill `0 <= m <= n*(n-1)//2`.

    samples : int, default=1
        Number of independent realizations of the `G(n, m)` random variable.

    random_state : int, RandomState instance or None, default=None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    method : {"auto", "tracking_selection", "reservoir_sampling", "pool"}, \
            default="auto"
        Determines which algorithm is used for sampling.
        See :class:`~sklearn.utils.random.sample_without_replacement` for
        more info.

    return_as : "adjacency_matrix" or "edge_array", default="edge_array"
        Determines whether the output is an an array containing pairs of
        vertices or a list of coo_matrices as samples of the adjacency matrix.

    dtype : dtype or str, default="int64"
        An integer data type. The value `n*(n-1)//2` must representable.
        No checks are implemented.

    Returns
    -------
    int ndarray or list of int coo_matrix:
        If return_as == "adjacency_matrix", a list of length `samples`
        containing sparse matrices of shape `(n, n)` with zeros and ones is
        returned;
        If return_as == "edge_array", an array of shape `(2, m, samples)`
        containing vertex indices in `[0, ..., n-1]` is returned.

    Examples
    --------
    >>> from sklearn.utils.random import sample_erdos_renyi_gnm
    >>> sample_erdos_renyi_gnm(5, 4, 3, random_state=42)
    array([[[4, 1, 4],
            [2, 2, 2],
            [3, 4, 1],
            [1, 3, 4]],
           [[2, 0, 3],
            [0, 0, 1],
            [2, 2, 0],
            [0, 2, 0]]])

    References
    ----------
    Erdős, P., Rényi, A. (1959).
    `On Random Graphs. <https://www.renyi.hu/~p_erdos/1959-11.pdf>`_
    Publicationes Mathematicae.
    """
    nedges = n*(n-1)//2

    # each sample has m edges, this array will eventually store their vertices
    if return_as == "adjacency_matrix":
        # leave room for conjugate elements
        er = np.empty([2, 2*m, samples], dtype=dtype)
    else:
        er = np.empty([2, m, samples], dtype=dtype)

    # make sure random state is valid
    random_state = check_random_state(random_state)

    # sample edge indices
    er[1, :m, :] = np.array(
        [sample_without_replacement(nedges, m, method=method,
                                    random_state=random_state)
            for _ in range(samples)], dtype=dtype).T

    # solves (t+1)*t//2 = er1+1 for t to get the left vertex index
    er[0, :m, :] = np.ceil(np.sqrt((er[1, :m, :]+1.)*2.+.25)-.5).astype(dtype)

    # computes er1 = er0-t*(t-1)//2 for the right vertex index
    er[1, :m, :] -= er[0, :m, :]*(er[0, :m, :]-1)//2

    if return_as == "adjacency_matrix":
        # set conjugate elements
        er[0, m:, :] = er[1, :m, :]
        er[1, m:, :] = er[0, :m, :]
        # set matrix values
        vals = np.ones(2*m, dtype=dtype)
        out = []
        # build a sparse matrix for each sample
        for k in range(samples):
            A = sp.coo_matrix(
                (vals, (er[0, :, k], er[1, :, k])), shape=(n, n), dtype=dtype)
            out.append(A)
        return out
    else:
        return er
