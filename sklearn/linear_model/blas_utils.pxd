# generic-type BLAS subroutines for l2 norm
from cython cimport floating

ctypedef floating (*BLAS_NRM2)(int n,
                               floating *w,
                               int incw) nogil

# generic-type BLAS subroutines for inplace scaling: w *= alpha
ctypedef void (*BLAS_SCAL)(int n,
                           floating alpha,
                           floating *w,
                           int incw) nogil

cdef extern from "cblas.h":
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    double dnrm2 "cblas_dnrm2"(int N,
                               double *X,
                               int incX) nogil
    float snrm2 "cblas_snrm2"(int N,
                              float *X,
                              int incX) nogil

    void dscal "cblas_dscal"(int N,
                             double alpha,
                             double *X,
                             int incX) nogil
    void sscal "cblas_sscal"(int N,
                             float alpha,
                             float *X,
                             int incX) nogil


