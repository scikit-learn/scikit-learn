from libc.math cimport exp, log
import numpy as np
cimport numpy as np
cimport cython

np.import_array()

ctypedef np.float64_t dtype_t

cdef dtype_t _NINF = -np.inf

@cython.boundscheck(False)
cdef dtype_t _max(dtype_t[:] values):
    # find maximum value (builtin 'max' is unrolled for speed)
    cdef dtype_t value
    cdef dtype_t vmax = _NINF
    for i in range(values.shape[0]):
        value = values[i]
        if value > vmax:
            vmax = value
    return vmax

@cython.boundscheck(False)
cpdef dtype_t _logsum(dtype_t[:] X):
    cdef dtype_t vmax = _max(X)
    cdef dtype_t power_sum = 0

    for i in range(X.shape[0]):
        power_sum += exp(X[i]-vmax)

    return log(power_sum) + vmax


@cython.boundscheck(False)
def _forward(int n_observations, int n_components,
        np.ndarray[dtype_t, ndim=1] log_startprob,
        np.ndarray[dtype_t, ndim=2] log_transmat,
        np.ndarray[dtype_t, ndim=2] framelogprob,
        np.ndarray[dtype_t, ndim=2] fwdlattice):

    cdef int t, i, j
    cdef double logprob
    cdef np.ndarray[dtype_t, ndim = 1] work_buffer
    work_buffer = np.zeros(n_components)

    for i in range(n_components):
        fwdlattice[0, i] = log_startprob[i] + framelogprob[0, i]

    for t in range(1, n_observations):
        for j in range(n_components):
            for i in range(n_components):
                work_buffer[i] = fwdlattice[t - 1, i] + log_transmat[i, j]
            fwdlattice[t, j] = _logsum(work_buffer) + framelogprob[t, j]


@cython.boundscheck(False)
def _backward(int n_observations, int n_components,
        np.ndarray[dtype_t, ndim=1] log_startprob,
        np.ndarray[dtype_t, ndim=2] log_transmat,
        np.ndarray[dtype_t, ndim=2] framelogprob,
        np.ndarray[dtype_t, ndim=2] bwdlattice):

    cdef int t, i, j
    cdef double logprob
    cdef np.ndarray[dtype_t, ndim = 1] work_buffer
    work_buffer = np.zeros(n_components)

    for i in range(n_components):
        bwdlattice[n_observations - 1, i] = 0.0

    for t in range(n_observations - 2, -1, -1):
        for i in range(n_components):
            for j in range(n_components):
                work_buffer[j] = log_transmat[i, j] + framelogprob[t + 1, j] \
                    + bwdlattice[t + 1, j]
            bwdlattice[t, i] = _logsum(work_buffer)


@cython.boundscheck(False)
def _compute_lneta(int n_observations, int n_components,
        np.ndarray[dtype_t, ndim=2] fwdlattice,
        np.ndarray[dtype_t, ndim=2] log_transmat,
        np.ndarray[dtype_t, ndim=2] bwdlattice,
        np.ndarray[dtype_t, ndim=2] framelogprob,
        double logprob,
        np.ndarray[dtype_t, ndim=3] lneta):

    cdef int i, j, t
    for t in range(n_observations - 1):
        for i in range(n_components):
            for j in range(n_components):
                lneta[t, i, j] = fwdlattice[t, i] + log_transmat[i, j] \
                    + framelogprob[t + 1, j] + bwdlattice[t + 1, j] - logprob


@cython.boundscheck(False)
def _viterbi(int n_observations, int n_components,
        np.ndarray[dtype_t, ndim=1] log_startprob,
        np.ndarray[dtype_t, ndim=2] log_transmat,
        np.ndarray[dtype_t, ndim=2] framelogprob):

    cdef int t, max_pos
    cdef np.ndarray[dtype_t, ndim = 2] viterbi_lattice
    cdef np.ndarray[np.int_t, ndim = 1] state_sequence
    cdef dtype_t logprob
    cdef np.ndarray[dtype_t, ndim = 2] work_buffer

    # Initialization
    state_sequence = np.empty(n_observations, dtype=np.int)
    viterbi_lattice = np.zeros((n_observations, n_components))
    viterbi_lattice[0] = log_startprob + framelogprob[0]

    # Induction
    for t in range(1, n_observations):
        work_buffer = viterbi_lattice[t-1] + log_transmat.T
        viterbi_lattice[t] = np.max(work_buffer, axis=1) + framelogprob[t]

    # Observation traceback
    max_pos = np.argmax(viterbi_lattice[n_observations - 1, :])
    state_sequence[n_observations - 1] = max_pos
    logprob = viterbi_lattice[n_observations - 1, max_pos]

    for t in range(n_observations - 2, -1, -1):
        max_pos = np.argmax(viterbi_lattice[t, :] \
                + log_transmat[:, state_sequence[t + 1]])
        state_sequence[t] = max_pos

    return state_sequence, logprob
