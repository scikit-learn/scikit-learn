cimport cython
import numpy as np
cimport numpy as np
from libc cimport math


cdef float EPSILON = 1e-16
cdef float PERPLEXITY_TOLERANCE = 1e-5


@cython.boundscheck(False)
cpdef np.ndarray[np.float_t, ndim=2] _binary_search_perplexity(
        np.ndarray[np.float_t, ndim=2] dist, float desired_perplexity,
        int verbose):
    """Binary search for sigmas of conditional Gaussians.

    We are looking for sigma = sqrt(1/beta) so that the perplexity
    of the conditional distribution p_i|j matches approximately the
    desired value. p_i|j are Gaussian distributed.

    The idea has been proposed in "Stochastic Neighbor Embedding"
    Geoffrey Hinton and Sam Roweis, 2003.

    Parameters
    ----------
    dist : array-like, shape (n_samples, n_samples)
        Distances between training samples.

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
    cdef int n_steps = 100

    cdef int n_samples = dist.shape[0]
    cdef np.ndarray[np.float_t, ndim=2] P = np.ndarray((n_samples, n_samples),
                                                       dtype=np.float)
    # Precisions of conditional Gaussian distrubutions
    cdef float beta
    cdef float beta_min
    cdef float beta_max
    cdef float beta_sum = 0.0
    # Now we go to log scale
    cdef float desired_entropy = math.log(desired_perplexity)
    cdef float entropy_diff

    cdef float entropy
    cdef float sum_Pi
    cdef float sum_disti_Pi
    cdef int i
    cdef int j

    for i in range(n_samples):
        beta_min = -math.INFINITY
        beta_max = math.INFINITY
        beta = 1.0

        # Binary search of precision for i-th conditional distribution
        for _ in range(n_steps):
            # Compute current entropy and corresponding probabilities
            for j in range(n_samples):
                P[i, j] = math.exp(-dist[i, j] * beta)
            P[i, i] = 0.0
            sum_Pi = 0.0
            for j in range(n_samples):
                sum_Pi += P[i, j]
            if sum_Pi == 0.0:
                sum_Pi = EPSILON
            sum_disti_Pi = 0.0
            for j in range(n_samples):
                P[i, j] /= sum_Pi
                sum_disti_Pi += dist[i, j] * P[i, j]
            entropy = math.log(sum_Pi) + beta * sum_disti_Pi
            entropy_diff = entropy - desired_entropy

            if math.fabs(entropy_diff) <= PERPLEXITY_TOLERANCE:
                break

            if entropy_diff > 0.0:
                beta_min = beta
                if beta_max == math.INFINITY:
                    beta *= 2.0
                else:
                    beta = (beta + beta_max) / 2.0
            else:
                beta_max = beta;
                if beta_min == -math.INFINITY:
                    beta /= 2.0
                else:
                    beta = (beta + beta_min) / 2.0

        beta_sum += beta

        if verbose and ((i + 1) % 1000 == 0 or i + 1 == n_samples):
            print("[t-SNE] Computed conditional probabilities for sample "
                  "%d / %d" % (i + 1, n_samples))

    if verbose:
        print("[t-SNE] Mean sigma: %f" % np.mean(np.sqrt(n_samples /
                                                         beta_sum)))
    return P
