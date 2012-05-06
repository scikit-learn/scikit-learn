import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h":
    double exp(double)
    double log(double)

ctypedef np.float64_t dtype_t


@cython.boundscheck(False)
def _logsum(int N, np.ndarray[dtype_t, ndim=1] X):
    cdef int i
    cdef double maxv, Xsum
    Xsum = 0.0
    maxv = X.max()
    for i in xrange(N):
        Xsum += exp(X[i] - maxv)
    return log(Xsum) + maxv


@cython.boundscheck(False)
def _forward(int n_observations, int n_components, \
        np.ndarray[dtype_t, ndim=1] log_startprob, \
        np.ndarray[dtype_t, ndim=2] log_transmat, \
        np.ndarray[dtype_t, ndim=2] framelogprob, \
        np.ndarray[dtype_t, ndim=2] fwdlattice):

    cdef int t, i, j
    cdef double logprob
    cdef np.ndarray[dtype_t, ndim = 1] work_buffer
    work_buffer = np.zeros(n_components)

    for i in xrange(n_components):
        fwdlattice[0, i] = log_startprob[i] + framelogprob[0, i]

    for t in xrange(1, n_observations):
        for j in xrange(n_components):
            for i in xrange(n_components):
                work_buffer[i] = fwdlattice[t - 1, i] + log_transmat[i, j]
            fwdlattice[t, j] = _logsum(n_components, work_buffer) \
                + framelogprob[t, j]


@cython.boundscheck(False)
def _backward(int n_observations, int n_components, \
        np.ndarray[dtype_t, ndim=1] log_startprob, \
        np.ndarray[dtype_t, ndim=2] log_transmat, \
        np.ndarray[dtype_t, ndim=2] framelogprob, \
        np.ndarray[dtype_t, ndim=2] bwdlattice):

    cdef int t, i, j
    cdef double logprob
    cdef np.ndarray[dtype_t, ndim = 1] work_buffer
    work_buffer = np.zeros(n_components)

    for i in xrange(n_components):
        bwdlattice[n_observations - 1, i] = 0.0

    for t in xrange(n_observations - 2, -1, -1):
        for i in xrange(n_components):
            for j in xrange(n_components):
                work_buffer[j] = log_transmat[i, j] + framelogprob[t + 1, j] \
                    + bwdlattice[t + 1, j]
            bwdlattice[t, i] = _logsum(n_components, work_buffer)


@cython.boundscheck(False)
def _compute_lneta(int n_observations, int n_components, \
        np.ndarray[dtype_t, ndim=2] fwdlattice, \
        np.ndarray[dtype_t, ndim=2] log_transmat, \
        np.ndarray[dtype_t, ndim=2] bwdlattice, \
        np.ndarray[dtype_t, ndim=2] framelogprob, \
        double logprob, \
        np.ndarray[dtype_t, ndim=3] lneta):

    cdef int i, j, t
    for t in xrange(n_observations - 1):
        for i in xrange(n_components):
            for j in xrange(n_components):
                lneta[t, i, j] = fwdlattice[t, i] + log_transmat[i, j] \
                    + framelogprob[t + 1, j] + bwdlattice[t + 1, j] - logprob


@cython.boundscheck(False)
def _viterbi(int n_observations, int n_components, \
        np.ndarray[dtype_t, ndim=1] log_startprob, \
        np.ndarray[dtype_t, ndim=2] log_transmat, \
        np.ndarray[dtype_t, ndim=2] framelogprob):

    cdef int i, j, t, max_pos
    cdef np.ndarray[dtype_t, ndim = 2] viterbi_lattice
    cdef np.ndarray[np.int_t, ndim = 1] state_sequence
    cdef double logprob
    cdef np.ndarray[dtype_t, ndim = 1] work_buffer

    # Initialize state_sequenceation
    state_sequence = np.zeros(n_observations, dtype=np.int)
    work_buffer = np.zeros(n_components)
    viterbi_lattice = np.zeros((n_observations, n_components))

    # viterbi_lattice[0,:] = log_startprob[:] + framelogprob[0,:]
    for i in xrange(n_components):
        viterbi_lattice[0, i] = log_startprob[i] + framelogprob[0, i]

    # Induction
    for t in xrange(1, n_observations):
        for j in xrange(n_components):
            work_buffer[:] = viterbi_lattice[t - 1, :] + log_transmat[:, j]
            viterbi_lattice[t, j] = np.max(work_buffer[:]) + framelogprob[t, j]

    # observation traceback
    max_pos = np.argmax(viterbi_lattice[n_observations - 1, :])
    state_sequence[n_observations - 1] = max_pos
    logprob = viterbi_lattice[n_observations - 1, max_pos]
    
    for t in xrange(n_observations - 2, -1, -1):
        max_pos = np.argmax(viterbi_lattice[t, :] \
                + log_transmat[:, state_sequence[t + 1]])
        state_sequence[t] = max_pos

    return state_sequence, logprob
