# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# Authors: Robert Layton <robertlayton@gmail.com>
#           Corey Lynch <coreylynch9@gmail.com>
# Licence: BSD 3 clause

from libc.math cimport exp
from scipy.special import gammaln
import numpy as np
cimport numpy as np
cimport cython
from sklearn.utils.lgamma cimport lgamma

np.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
def expected_mutual_information(contingency, int n_samples):
    """Calculate the expected mutual information for two labelings."""
    cdef int R, C
    cdef float N, gln_N, emi, term2, term3, gln
    cdef np.ndarray[double] gln_a, gln_b, gln_Na, gln_Nb, gln_nij, log_Nnij
    cdef np.ndarray[double] nijs, term1
    cdef np.ndarray[double, ndim=2] log_ab_outer
    cdef np.ndarray[np.int32_t] a, b
    #cdef np.ndarray[int, ndim=2] start, end
    R, C = contingency.shape
    N = float(n_samples)
    a = np.sum(contingency, axis=1).astype(np.int32)
    b = np.sum(contingency, axis=0).astype(np.int32)
    # There are three major terms to the EMI equation, which are multiplied to
    # and then summed over varying nij values.
    # While nijs[0] will never be used, having it simplifies the indexing.
    nijs = np.arange(0, max(np.max(a), np.max(b)) + 1, dtype='float')
    nijs[0] = 1  # Stops divide by zero warnings. As its not used, no issue.
    # term1 is nij / N
    term1 = nijs / N
    # term2 is log((N*nij) / (a * b)) == log(N * nij) - log(a * b)
    # term2 uses the outer product
    log_ab_outer = np.log(a)[:, np.newaxis] + np.log(b)
    # term2 uses N * nij
    log_Nnij = np.log(N * nijs)
    # term3 is large, and involved many factorials. Calculate these in log
    # space to stop overflows.
    gln_a = gammaln(a + 1)
    gln_b = gammaln(b + 1)
    gln_Na = gammaln(N - a + 1)
    gln_Nb = gammaln(N - b + 1)
    gln_N = gammaln(N + 1)
    gln_nij = gammaln(nijs + 1)
    # start and end values for nij terms for each summation.
    start = np.array([[v - N + w for w in b] for v in a], dtype='int')
    start = np.maximum(start, 1)
    end = np.minimum(np.resize(a, (C, R)).T, np.resize(b, (R, C))) + 1
    # emi itself is a summation over the various values.
    emi = 0
    cdef Py_ssize_t i, j, nij
    for i in range(R):
        for j in range(C):
            for nij in range(start[i,j], end[i,j]):
                term2 = log_Nnij[nij] - log_ab_outer[i,j]
                # Numerators are positive, denominators are negative.
                gln = (gln_a[i] + gln_b[j] + gln_Na[i] + gln_Nb[j]
                     - gln_N - gln_nij[nij] - lgamma(a[i] - nij + 1)
                     - lgamma(b[j] - nij + 1)
                     - lgamma(N - a[i] - b[j] + nij + 1))
                term3 = exp(gln)
                emi += (term1[nij] * term2 * term3)
    return emi
