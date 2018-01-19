# Synopsis: Unified consistent BLAS wrappers
# Author: Elvis Dohmatob <gmdopp@gmail.com>
#
# Notes: A modern alternative would be to use scipy.linalg.cython_blas

from cython cimport floating

cdef extern from "cblas.h" nogil:
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    # LEVEL 1: in-place scaling
    void sscal "cblas_sscal"(int N, float alpha, float *X, int incX)
    void dscal "cblas_dscal"(int N, double alpha, double *X, int incX)

    # LEVEL 1: sum of absolute values
    float sasum "cblas_sasum"(int N, float *X, int incX)
    double dasum "cblas_dasum"(int N, double *X, int incX)

    # LEVEL 1: copy from one array to another
    void scopy "cblas_scopy"(int N, float *X, int incX, float *Y, int incY)
    void dcopy "cblas_dcopy"(int N, double *X, int incX, double *Y, int incY)

    # LEVEL 1: vector update: Y = Y + alpha * X
    void saxpy "cblas_saxpy"(int N, float alpha, float *X, int incX, float *Y,
                             int incY)
    void daxpy "cblas_daxpy"(int N, double alpha, double *X, int incX, double *Y,
                             int incY)

    # LEVEL 2: L2 norm
    float snrm2 "cblas_snrm2"(int N, float *X, int incX)
    double dnrm2 "cblas_dnrm2"(int N, double *X, int incX)

    # LEVEL 2: dot products
    float sdot "cblas_sdot"(int N, float *X, int incX, float *Y, int incY)
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY)

    # LEVEL 3: rank-1 updates
    void sger "cblas_sger"(CBLAS_ORDER Order, int M, int N, float alpha,
                           float *X, int incX, float *Y, int incY, float *A,
                           int lda)
    void dger "cblas_dger"(CBLAS_ORDER Order, int M, int N, double alpha,
                           double *X, int incX, double *Y, int incY, double *A,
                           int lda)

### Wrap BLAS subroutines into a consistent API

# copy data from one array (X) to another (Y)
cdef void fused_copy (int N, floating *X, int incX, floating *Y,
                      int incY) nogil

# multiply a array (X) with a scalar (alpha) in-place
cdef void fused_scal(int N, floating alpha, floating *X, int incX) nogil

# compute some of absolute values of an array (X)
cdef floating fused_asum(int N, floating *X, int incX) nogil

# vector update: Y = Y + alpha * X
cdef void fused_axpy(int N, floating alpha, floating *X, int incX,
                     floating *Y, int incY) nogil

# computes dot inner product of X with Y
cdef floating fused_dot(int N, floating *X, int incX, floating *Y,
                           int incY) nogil

# rank-1 update: A = A + alpha * outer(X, Y)
cdef void fused_ger(CBLAS_ORDER Order, int M, int N, floating alpha,
                    floating *X, int incX, floating *Y, int incY,
                    floating *A, int lda) nogil

# norm of array: normX = ||X||_2
cdef floating fused_nrm2(int N, floating *X, int incX) nogil
