# cython: language_level=3, boundscheck=False

from libc cimport math
cimport cython
import numpy as np
cimport numpy as np
from libc.stdio cimport printf
cdef extern from "numpy/npy_math.h":
    float NPY_INFINITY


cdef float EPSILON_DBL = 1e-8
cdef float PERPLEXITY_TOLERANCE = 1e-5

cpdef np.ndarray[np.float32_t, ndim=2] _binary_search_perplexity(
        np.ndarray[np.float32_t, ndim=2] affinities,
        np.ndarray[np.int64_t, ndim=2] neighbors,
        float desired_perplexity,
        int verbose):
    """Binary search for sigmas of conditional Gaussians.

    This approximation reduces the computational complexity from O(N^2) to
    O(uN). See the exact method '_binary_search_perplexity' for more details.

    Parameters
    ----------
    affinities : array-like, shape (n_samples, k)
        Distances between training samples and its k nearest neighbors.

    neighbors : array-like, shape (n_samples, k) or None
        Each row contains the indices to the k nearest neigbors. If this
        array is None, then the perplexity is estimated over all data
        not just the nearest neighbors.

    desired_perplexity : float
        Desired perplexity (2^entropy) of the conditional Gaussians.

    verbose : int
        Verbosity level.

    Returns
    -------
    P : array, shape (n_samples, n_samples)
        Probabilities of conditional Gaussian distributions p_i|j.
    """
    # Maximum number of binary search steps
    cdef long n_steps = 100

    cdef long n_samples = affinities.shape[0]
    # Precisions of conditional Gaussian distributions
    cdef float beta
    cdef float beta_min
    cdef float beta_max
    cdef float beta_sum = 0.0

    # Use log scale
    cdef float desired_entropy = math.log(desired_perplexity)
    cdef float entropy_diff

    cdef float entropy
    cdef float sum_Pi
    cdef float sum_disti_Pi
    cdef long i, j, k, l
    cdef long n_neighbors = n_samples
    cdef int using_neighbors = neighbors is not None

    if using_neighbors:
        n_neighbors = neighbors.shape[1]

    # This array is later used as a 32bit array. It has multiple intermediate
    # floating point additions that benefit from the extra precision
    cdef np.ndarray[np.float64_t, ndim=2] P = np.zeros(
        (n_samples, n_neighbors), dtype=np.float64)

    for i in range(n_samples):
        beta_min = -NPY_INFINITY
        beta_max = NPY_INFINITY
        beta = 1.0

        # Binary search of precision for i-th conditional distribution
        for l in range(n_steps):
            # Compute current entropy and corresponding probabilities
            # computed just over the nearest neighbors or over all data
            # if we're not using neighbors
            sum_Pi = 0.0
            for j in range(n_neighbors):
                if j != i or using_neighbors:
                    P[i, j] = math.exp(-affinities[i, j] * beta)
                    sum_Pi += P[i, j]

            if sum_Pi == 0.0:
                sum_Pi = EPSILON_DBL
            sum_disti_Pi = 0.0

            for j in range(n_neighbors):
                P[i, j] /= sum_Pi
                sum_disti_Pi += affinities[i, j] * P[i, j]

            entropy = math.log(sum_Pi) + beta * sum_disti_Pi
            entropy_diff = entropy - desired_entropy

            if math.fabs(entropy_diff) <= PERPLEXITY_TOLERANCE:
                break

            if entropy_diff > 0.0:
                beta_min = beta
                if beta_max == NPY_INFINITY:
                    beta *= 2.0
                else:
                    beta = (beta + beta_max) / 2.0
            else:
                beta_max = beta
                if beta_min == -NPY_INFINITY:
                    beta /= 2.0
                else:
                    beta = (beta + beta_min) / 2.0

        beta_sum += beta

        if verbose and ((i + 1) % 1000 == 0 or i + 1 == n_samples):
            print("[t-SNE] Computed conditional probabilities for sample "
                  "%d / %d" % (i + 1, n_samples))

    if verbose:
        print("[t-SNE] Mean sigma: %f"
              % np.mean(math.sqrt(n_samples / beta_sum)))
    return P
