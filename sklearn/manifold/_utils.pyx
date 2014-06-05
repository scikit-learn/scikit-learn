from libc cimport math
cimport cython
import numpy as np
cimport numpy as np
cdef extern from "numpy/npy_math.h":
    float NPY_INFINITY


cdef double EPSILON_DBL = 1e-7
cdef double PERPLEXITY_TOLERANCE = 1e-5


@cython.boundscheck(False)
cpdef np.ndarray[np.float_t, ndim=2] _binary_search_perplexity(
        np.ndarray[np.float_t, ndim=2] affinities, double desired_perplexity,
        int verbose):
    """Binary search for sigmas of conditional Gaussians.

    We are looking for sigma = sqrt(1/beta) so that the perplexity
    of the conditional distribution p_i|j matches approximately the
    desired value. p_i|j are Gaussian distributed.

    The idea has been proposed in "Stochastic Neighbor Embedding"
    Geoffrey Hinton and Sam Roweis, 2003.

    Parameters
    ----------
    affinities : array-like, shape (n_samples, n_samples)
        Distances between training samples.

    desired_perplexity : double
        Desired perplexity (2^entropy) of the conditional Gaussians.

    verbose : int
        Verbosity level.

    Returns
    -------
    P : array, shape (n_samples, n_samples)
        Probabilities of conditional Gaussian distributions p_i|j.
    """
    # Maximum number of binary search steps
    cdef int n_steps = 100

    cdef int n_samples = affinities.shape[0]
    cdef np.ndarray[np.float_t, ndim=2] P = np.ndarray((n_samples, n_samples),
                                                       dtype=np.double)
    # Precisions of conditional Gaussian distrubutions
    cdef double beta
    cdef double beta_min
    cdef double beta_max
    cdef double beta_sum = 0.0
    # Now we go to log scale
    cdef double desired_entropy = math.log(desired_perplexity)
    cdef double entropy_diff

    cdef double entropy
    cdef double sum_Pi
    cdef double sum_disti_Pi
    cdef int i
    cdef int j

    for i in range(n_samples):
        beta_min = -NPY_INFINITY
        beta_max = NPY_INFINITY
        beta = 1.0

        # Binary search of precision for i-th conditional distribution
        for _ in range(n_steps):
            # Compute current entropy and corresponding probabilities
            for j in range(n_samples):
                P[i, j] = math.exp(-affinities[i, j] * beta)
            P[i, i] = 0.0
            sum_Pi = 0.0
            for j in range(n_samples):
                sum_Pi += P[i, j]
            if sum_Pi == 0.0:
                sum_Pi = EPSILON_DBL
            sum_disti_Pi = 0.0
            for j in range(n_samples):
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
