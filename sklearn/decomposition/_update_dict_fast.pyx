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


np.import_array()


cdef extern from "numpy/arrayobject.h":
    bint PyArray_IS_C_CONTIGUOUS(ndarray)
    bint PyArray_IS_F_CONTIGUOUS(ndarray)


cdef extern from "cblas.h":
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y,
                             int incY) nogil
    float sdot "cblas_sdot"(int N, float *X, int incX, float *Y,
                            int incY) nogil
    void dgemm "cblas_dgemm"(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                             CBLAS_TRANSPOSE TransB, int M, int N,
                             int K, double alpha, double *A,
                             int lda, double *B, int ldb,
                             double beta, double *C, int ldc) nogil
    void sgemm "cblas_sgemm"(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                             CBLAS_TRANSPOSE TransB, int M, int N,
                             int K, float alpha, float *A,
                             int lda, float *B, int ldb,
                             float beta, float *C, int ldc) nogil
    void dgemv "cblas_dgemv"(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                             int M, int N, double alpha, double *A, int lda,
                             double *X, int incX, double beta,
                             double *Y, int incY) nogil
    void sgemv "cblas_sgemv"(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                             int M, int N, float alpha, float *A, int lda,
                             float *X, int incX, float beta,
                             float *Y, int incY) nogil
    void dger "cblas_dger"(CBLAS_ORDER Order, int M, int N, double alpha,
                           double *X, int incX, double *Y, int incY,
                           double *A, int lda) nogil
    void sger "cblas_sger"(CBLAS_ORDER Order, int M, int N, float alpha,
                           float *X, int incX, float *Y, int incY,
                           float *A, int lda) nogil
    void dscal "cblas_dscal"(int N, double alpha, double *X, int incX) nogil
    void sscal "cblas_sscal"(int N, float alpha, float *X, int incX) nogil


cdef inline floating dot(int N,
                         floating* X, int incX,
                         floating* Y, int incY) nogil:
    if floating is float:
        return sdot(N, X, incX, Y, incY)
    else:
        return ddot(N, X, incX, Y, incY)


cdef inline void gemm(CBLAS_ORDER Order,
                      CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
                      int M, int N, int K,
                      floating alpha, floating* A, int lda,
                      floating* B, int ldb,
                      floating beta, floating* C, int ldc) nogil:
    if floating is float:
        sgemm(Order, TransA, TransB,
              M, N, K,
              alpha, A, lda,
              B, ldb,
              beta, C, ldc)
    else:
        dgemm(Order, TransA, TransB,
              M, N, K,
              alpha, A, lda,
              B, ldb,
              beta, C, ldc)


cdef inline void gemv(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                      int M, int N,
                      floating alpha, floating *A, int lda,
                      floating *X, int incX,
                      floating beta, floating *Y, int incY) nogil:
    if floating is float:
        sgemv(Order, TransA,
              M, N,
              alpha, A, lda,
              X, incX,
              beta, Y, incY)
    else:
        dgemv(Order, TransA,
              M, N,
              alpha, A, lda,
              X, incX,
              beta, Y, incY)


cdef inline void ger(CBLAS_ORDER Order,
                     int M, int N,
                     floating alpha, floating *X, int incX,
                     floating *Y, int incY,
                     floating *A, int lda) nogil:
    if floating is float:
        sger(Order,
             M, N,
             alpha, X, incX,
             Y, incY,
             A, lda)
    else:
        dger(Order,
             M, N,
             alpha, X, incX,
             Y, incY,
             A, lda)


cdef inline void scal(int N, floating alpha, floating* X, int incX) nogil:
    if floating is float:
        sscal(N, alpha, X, incX)
    else:
        dscal(N, alpha, X, incX)


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
cdef inline void randn_atom_k(uint32_t* seed,
                              npy_intp k,
                              floating[:, :] a) nogil:
    cdef npy_intp i, r, n
    cdef floating tmp

    n = a.shape[0]
    r = n % 2

    for i in range(0, n - r, 2):
        our_randn(seed, &a[i, k], &a[i + 1, k])
    if r == 1:
        our_randn(seed, &a[n - 1, k], &tmp)


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline void clip_negative_k(npy_intp k, floating[:, :] a) nogil:
    cdef npy_intp i
    for i in range(a.shape[0]):
        a[i, k] = fmax(a[i, k], 0.0)


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

    # Get bounds
    cdef npy_intp n_features
    cdef npy_intp n_samples
    cdef npy_intp n_components

    # Indices to iterate over
    cdef npy_intp k

    # For holding the kth atom
    cdef floating atom_norm2

    # Verbose message
    cdef char* msg

    # Random number seed for use in C
    cdef uint32_t rand_r_state
    cdef uint32_t* rand_r_state_ptr

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
        # Assign bounds
        n_features = Y.shape[0]
        n_components = code_C.shape[0]
        n_samples = Y.shape[1]

        # Denote whether an IO error occurred
        ioerr = False

        # Determine verbose message
        if verbose == 0:
            msg = NULL
        elif verbose == 1:
            msg = b"+"
        else:
            msg = b"Adding new random atom"

        # Get pointer to random state
        rand_r_state_ptr = &rand_r_state

        # R <- -1.0 * U * V^T + 1.0 * Y
        gemm(CblasColMajor, CblasNoTrans, CblasTrans,
             n_features, n_samples, n_components,
             -1.0, &dictionary_F[0, 0], n_features,
             &code_C[0, 0], n_samples,
             1.0, &R[0, 0], n_features)

        for k in range(n_components):
            # R <- 1.0 * U_k * V_k^T + R
            ger(CblasColMajor, n_features, n_samples,
                1.0, &dictionary_F[0, k], 1,
                &code_C[k, 0], 1,
                &R[0, 0], n_features)

            # U_k <- 1.0 * R * V_k^T
            gemv(CblasColMajor, CblasNoTrans,
                 n_features, n_samples,
                 1.0, &R[0, 0], n_features,
                 &code_C[k, 0], 1,
                 0.0, &dictionary_F[0, k], 1)

            # Clip negative values
            if positive:
                clip_negative_k(k, dictionary_F)

            # Scale k'th atom
            # U_k * U_k
            atom_norm2 = dot(n_features,
                             &dictionary_F[0, k], 1,
                             &dictionary_F[0, k], 1)

            # Generate random atom to replace inconsequential one
            if atom_norm2 < 1e-20:
                # Handle verbose mode
                if msg is not NULL:
                    if puts(msg) == EOF or fflush(stdout) == EOF:
                        ioerr = True
                        break

                # Seed random atom
                randn_atom_k(rand_r_state_ptr, k, dictionary_F)

                # Clip negative values
                if positive:
                    clip_negative_k(k, dictionary_F)

                # Setting corresponding coefs to 0
                scal(n_samples, 0.0, &code_C[k, 0], 1)

                # Compute new norm
                # U_k * U_k
                atom_norm2 = dot(n_features,
                                 &dictionary_F[0, k], 1,
                                 &dictionary_F[0, k], 1)

                # Normalize atom
                scal(n_features,
                     1.0 / sqrt(atom_norm2),
                     &dictionary_F[0, k], 1)
            else:
                # Normalize atom
                scal(n_features,
                     1.0 / sqrt(atom_norm2),
                     &dictionary_F[0, k], 1)

                # R <- -1.0 * U_k * V_k^T + R
                ger(CblasColMajor, n_features, n_samples,
                    -1.0, &dictionary_F[0, k], 1,
                    &code_C[k, 0], 1,
                    &R[0, 0], n_features)

        # Compute sum of squared residuals
        if not ioerr and return_r2:
            R2 = dot(n_features * n_samples,
                     &R[0, 0], 1,
                     &R[0, 0], 1)

    # Raise if verbose printing failed
    if ioerr:
        raise IOError("Failed to print out state.")

    # Optionally return residuals
    if return_r2:
        return dictionary_arr, R2
    else:
        return dictionary_arr
