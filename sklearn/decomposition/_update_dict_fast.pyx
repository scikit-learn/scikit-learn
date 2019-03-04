cimport libc
cimport libc.math
cimport libc.stdio
from libc.math cimport fmax, sqrt, log as ln
from libc.stdio cimport EOF, fflush, puts, stdout

cimport cython
from cython cimport floating

cimport numpy as np
import numpy as np
from numpy cimport npy_intp, uint32_t

from sklearn.utils import check_random_state

from ..utils._cython_blas cimport _dot, _scal, _ger, _gemv, _gemm
from ..utils._cython_blas cimport RowMajor, ColMajor, Trans, NoTrans


np.import_array()


cdef struct bint_false_type:
    char empty

cdef struct bint_true_type:
    char empty


cdef extern from "numpy/arrayobject.h":
    bint PyArray_IS_C_CONTIGUOUS(ndarray)
    bint PyArray_IS_F_CONTIGUOUS(ndarray)


cdef inline np.ndarray ensure_c(np.ndarray arr, bint copy):
    if copy or not PyArray_IS_C_CONTIGUOUS(arr):
        return np.PyArray_NewCopy(arr, np.NPY_CORDER)
    else:
        return arr


cdef inline np.ndarray ensure_fortran(np.ndarray arr, bint copy):
    if copy or not PyArray_IS_F_CONTIGUOUS(arr):
        return np.PyArray_NewCopy(arr, np.NPY_FORTRANORDER)
    else:
        return arr


cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline uint32_t our_rand_r(uint32_t* seed) nogil:
    seed[0] ^= <uint32_t>(seed[0] << 13)
    seed[0] ^= <uint32_t>(seed[0] >> 17)
    seed[0] ^= <uint32_t>(seed[0] << 5)

    return seed[0] % (<uint32_t>RAND_R_MAX + 1)


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline void our_rand_float(uint32_t* seed, floating* out) nogil:
    cdef int i, n
    cdef floating x

    if floating is float:
        n = 1
    else:
        n = 2

    x = 0.0
    for i in range(n):
        x += our_rand_r(seed)
        x /= (<uint32_t>RAND_R_MAX + 1)

    out[0] = x


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef void our_randn(uint32_t* seed,
                    floating* out1,
                    floating* out2) nogil:
    cdef int i
    cdef floating s
    cdef floating x[2]

    while True:
        s = 0.0
        for i in range(2):
            our_rand_float(seed, &x[i])

            x[i] *= 2.0
            x[i] -= 1.0

            s += x[i] ** 2.0

        if 0.0 < s < 1.0:
            break

    s = sqrt((-2.0 * ln(s)) / s)

    x[0] *= s
    x[1] *= s

    out1[0] = x[0]
    out2[0] = x[1]


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline void randn_atom(uint32_t* seed, floating* a, npy_intp n) nogil:
    cdef npy_intp i, r
    cdef floating tmp

    r = n % 2

    for i in range(0, n - r, 2):
        our_randn(seed, &a[i], &a[i + 1])
    if r == 1:
        our_randn(seed, &a[n - 1], &tmp)


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline void clip_negative(floating* a, npy_intp n) nogil:
    cdef npy_intp i
    for i in range(n):
        a[i] = fmax(a[i], 0.0)


cdef fused positive_bint_type:
    bint_false_type
    bint_true_type


cdef fused verbose_bint_type:
    bint_false_type
    bint_true_type


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef bint compute_dict(npy_intp n_features,
                       npy_intp n_samples,
                       npy_intp n_components,
                       floating[::, :] dictionary_F,
                       floating[::, :] R,
                       floating[:, ::] code_C,
                       const char* msg,
                       uint32_t* rand_r_state_ptr,
                       floating* R2_ptr,
                       bint return_r2,
                       positive_bint_type* positive,
                       verbose_bint_type* verbose) nogil:

    # Indices to iterate over
    cdef npy_intp k

    # Pointers to data for readability
    cdef floating* dictionary_F_ptr
    cdef floating* code_C_ptr
    cdef floating* R_ptr

    # Pointers to kth atom and kth sparse code
    cdef floating* kth_dictionary_F_ptr
    cdef floating* kth_code_C_ptr

    # For holding the kth atom
    cdef floating atom_norm2

    # Get pointer to memoryviews' data
    dictionary_F_ptr = &dictionary_F[0, 0]
    code_C_ptr = &code_C[0, 0]
    R_ptr = &R[0, 0]

    # R <- -1.0 * U * V^T + 1.0 * Y
    _gemm(ColMajor, NoTrans, Trans,
          n_features, n_samples, n_components,
          -1.0, dictionary_F_ptr, n_features,
          code_C_ptr, n_samples,
          1.0, R_ptr, n_features)

    for k in range(n_components):
        # Get pointers to current atom and code
        kth_dictionary_F_ptr = &dictionary_F[0, k]
        kth_code_C_ptr = &code_C[k, 0]

        # R <- 1.0 * U_k * V_k^T + R
        _ger(ColMajor, n_features, n_samples,
             1.0, kth_dictionary_F_ptr, 1,
             kth_code_C_ptr, 1,
             R_ptr, n_features)

        # U_k <- 1.0 * R * V_k^T
        _gemv(ColMajor, NoTrans,
              n_features, n_samples,
              1.0, R_ptr, n_features,
              kth_code_C_ptr, 1,
              0.0, kth_dictionary_F_ptr, 1)

        # Clip negative values
        if positive_bint_type is bint_true_type:
            clip_negative(kth_dictionary_F_ptr, n_features)

        # Scale k'th atom
        # U_k * U_k
        atom_norm2 = _dot(n_features,
                          kth_dictionary_F_ptr, 1,
                          kth_dictionary_F_ptr, 1)

        # Generate random atom to replace inconsequential one
        if atom_norm2 < 1e-20:
            # Handle verbose mode
            if verbose_bint_type is bint_true_type:
                if puts(msg) == EOF or fflush(stdout) == EOF:
                    return True

            # Seed random atom
            randn_atom(rand_r_state_ptr, kth_dictionary_F_ptr, n_features)

            # Clip negative values
            if positive_bint_type is bint_true_type:
                clip_negative(kth_dictionary_F_ptr, n_features)

            # Setting corresponding coefs to 0
            _scal(n_samples, 0.0, kth_code_C_ptr, 1)

            # Compute new norm
            # U_k * U_k
            atom_norm2 = _dot(n_features,
                              kth_dictionary_F_ptr, 1,
                              kth_dictionary_F_ptr, 1)

            # Normalize atom
            _scal(n_features, 1.0 / sqrt(atom_norm2), kth_dictionary_F_ptr, 1)
        else:
            # Normalize atom
            _scal(n_features, 1.0 / sqrt(atom_norm2), kth_dictionary_F_ptr, 1)

            # R <- -1.0 * U_k * V_k^T + R
            _ger(ColMajor, n_features, n_samples,
                 -1.0, kth_dictionary_F_ptr, 1,
                 kth_code_C_ptr, 1,
                 R_ptr, n_features)

    # Compute sum of squared residuals
    if return_r2:
        R2_ptr[0] = _dot(n_features * n_samples, R_ptr, 1, R_ptr, 1)

    return False


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline bint dispatch_compute_dict(npy_intp n_features,
                                       npy_intp n_samples,
                                       npy_intp n_components,
                                       floating[::, :] dictionary_F,
                                       floating[::, :] R,
                                       floating[:, ::] code_C,
                                       const char* msg,
                                       uint32_t* rand_r_state_ptr,
                                       floating* R2_ptr,
                                       bint positive,
                                       bint verbose,
                                       bint return_r2) nogil:

    # Positive/Verbose condition case
    cdef unsigned int case

    # Determine case case
    case = 2 * positive + verbose

    # Dispatch to specialization with right condition
    # Avoids repeated checks within a for-loop
    if case == 3:
        return compute_dict(n_features, n_samples, n_components,
                            dictionary_F, R, code_C,
                            msg, rand_r_state_ptr,
                            R2_ptr, return_r2,
                            <bint_true_type*>NULL,
                            <bint_true_type*>NULL)
    elif case == 2:
        return compute_dict(n_features, n_samples, n_components,
                            dictionary_F, R, code_C,
                            msg, rand_r_state_ptr,
                            R2_ptr, return_r2,
                            <bint_true_type*>NULL,
                            <bint_false_type*>NULL)
    elif case == 1:
        return compute_dict(n_features, n_samples, n_components,
                            dictionary_F, R, code_C,
                            msg, rand_r_state_ptr,
                            R2_ptr, return_r2,
                            <bint_false_type*>NULL,
                            <bint_true_type*>NULL)
    else:
        return compute_dict(n_features, n_samples, n_components,
                            dictionary_F, R, code_C,
                            msg, rand_r_state_ptr,
                            R2_ptr, return_r2,
                            <bint_false_type*>NULL,
                            <bint_false_type*>NULL)


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
def update_dict(np.ndarray[floating, ndim=2] dictionary not None,
                np.ndarray[floating, ndim=2] Y not None,
                np.ndarray[floating, ndim=2] code not None,
                unsigned int verbose=0, bint return_r2=False,
                random_state=None, bint positive=False):
    """Update the dense dictionary factor in place.

    Parameters
    ----------
    dictionary : array of shape (n_features, n_components)
        Value of the dictionary at the previous iteration.

    Y : array of shape (n_features, n_samples)
        Data matrix.

    code : array of shape (n_components, n_samples)
        Sparse coding of the data against which to optimize the dictionary.

    verbose:
        Degree of output the procedure will print.

    return_r2 : bool
        Whether to compute and return the residual sum of squares corresponding
        to the computed solution.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    positive : boolean, optional
        Whether to enforce positivity when finding the dictionary.

        .. versionadded:: 0.20

    Returns
    -------
    dictionary : array of shape (n_features, n_components)
        Updated dictionary.

    """

    # Array views/copies of the data
    cdef floating[::, :] dictionary_F
    cdef floating[:, ::] code_C
    cdef floating[::, :] R

    # Verbose message
    cdef char* msg

    # Random number seed for use in C
    cdef uint32_t rand_r_state

    # Results
    cdef bint ioerr
    cdef np.ndarray dictionary_arr
    cdef floating R2

    # Create random number seed to use
    random_state = check_random_state(random_state)
    rand_r_state = random_state.randint(0, RAND_R_MAX)

    # Initialize Contiguous Arrays
    dictionary_arr = ensure_fortran(dictionary, False)
    dictionary_F = dictionary_arr
    code_C = ensure_c(code, False)
    R = ensure_fortran(Y, True)

    with nogil:
        # Determine verbose message
        if verbose == 0:
            msg = NULL
        elif verbose == 1:
            msg = b"+"
        else:
            msg = b"Adding new random atom"

        # Compute updated dictionary
        ioerr = dispatch_compute_dict(Y.shape[0], Y.shape[1], code_C.shape[0],
                                      dictionary_F, R, code_C,
                                      msg, &rand_r_state, &R2,
                                      positive, verbose, return_r2)

    # Raise if verbose printing failed
    if ioerr:
        raise IOError("Failed to print out state.")

    # Optionally return residuals
    if return_r2:
        return dictionary_arr, R2
    else:
        return dictionary_arr
