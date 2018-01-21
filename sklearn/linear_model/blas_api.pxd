# Synopsis: Unified consistent BLAS wrappers
# Author: Elvis Dohmatob <gmdopp@gmail.com>
#
# Notes: A modern alternative would be to use scipy.linalg.cython_blas

from types cimport floating, complexing

cdef extern from "cblas.h" nogil:
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    # LEVEL 1
    void sscal "cblas_sscal"(int N,
                             float alpha,
                             float *X,
                             int incX)
    void dscal "cblas_dscal"(int N,
                             double alpha,
                             double *X,
                             int incX)
    void cscal "cblas_cscal"(int N,
                             void *alpha,
                             void *X,
                             int incX)
    void zscal "cblas_zscal"(int N,
                             void *alpha,
                             void *X,
                             int incX)

    float sasum "cblas_sasum"(int N,
                              float *X,
                              int incX)
    double dasum "cblas_dasum"(int N,
                               double *X,
                               int incX)

    void scopy "cblas_scopy"(int N,
                             float *X,
                             int incX,
                             float *Y,
                             int incY)
    void dcopy "cblas_dcopy"(int N,
                             double *X,
                             int incX,
                             double *Y,
                             int incY)
    void ccopy "cblas_ccopy"(int N,
                             void *X,
                             int incX,
                             void *Y,
                             int incY)
    void zcopy "cblas_zcopy"(int N,
                             void *X,
                             int incX,
                             void *Y,
                             int incY)

    void saxpy "cblas_saxpy"(int N,
                             float alpha,
                             float *X,
                             int incX,
                             float *Y,
                             int incY)
    void daxpy "cblas_daxpy"(int N,
                             double alpha,
                             double *X,
                             int incX,
                             double *Y,
                             int incY)
    void caxpy "cblas_caxpy"(int N,
                             void *alpha,
                             void *X,
                             int incX,
                             void *Y,
                             int incY)
    void zaxpy "cblas_zaxpy"(int N,
                             void *alpha,
                             void *X,
                             int incX,
                             void *Y,
                             int incY)



    # LEVEL 2
    float snrm2 "cblas_snrm2"(int N,
                              float *X,
                              int incX)
    double dnrm2 "cblas_dnrm2"(int N,
                               double *X,
                               int incX)
    double scnrm2 "cblas_scnrm2"(int N,
                                 void *X,
                                 int incX)
    double dznrm2 "cblas_dznrm2"(int N,
                                 void *X,
                                 int incX)

    float sdot "cblas_sdot"(int N,
                            float *X,
                            int incX,
                            float *Y,
                            int incY)
    double ddot "cblas_ddot"(int N,
                             double *X,
                             int incX,
                             double *Y,
                             int incY)
    float cdotu "cblas_cdotu_sub"(int N,
                                  void *X,
                                  int incX,
                                  void *Y,
                                  int incY,
                                  void *z)
    double zdotu "cblas_zdotu_sub"(int N,
                                   void *X,
                                   int incX,
                                   void *Y,
                                   int incY,
                                   void *z)
    float cdotc "cblas_cdotc_sub"(int N,
                                  void *X,
                                  int incX,
                                  void *Y,
                                  int incY,
                                  void *z)
    double zdotc "cblas_zdotc_sub"(int N,
                                   void *X,
                                   int incX,
                                   void *Y,
                                   int incY,
                                   void *z)

    # LEVEL 3
    void sger "cblas_sger"(CBLAS_ORDER Order,
                           int M,
                           int N,
                           float alpha,
                           float *X,
                           int incX,
                           float *Y,
                           int incY,
                           float *A,
                           int lda)
    void dger "cblas_dger"(CBLAS_ORDER Order,
                           int M,
                           int N,
                           double alpha,
                           double *X,
                           int incX,
                           double *Y,
                           int incY,
                           double *A,
                           int lda)
    void cgeru "cblas_cgeru"(CBLAS_ORDER Order,
                             int M,
                             int N,
                             void *alpha,
                             void *X,
                             int incX,
                             void *Y,
                             int incY,
                             void *A,
                             int lda)
    void zgeru "cblas_zgeru"(CBLAS_ORDER Order,
                             int M,
                             int N,
                             void *alpha,
                             void *X,
                             int incX,
                             void *Y,
                             int incY,
                             void *A,
                             int lda)

# Special treatment of trick BLAS APIs
cdef void fused_copy (int N,
                      complexing *X,
                      int incX,
                      complexing *Y,
                      int incY) nogil
cdef void fused_scal(int N,
                     complexing alpha,
                     complexing *X,
                     int incX) nogil
cdef double fused_asum(int N,
                       complexing *X,
                       int incX) nogil
cdef complexing fused_axpy(int N,
                           complexing alpha,
                           complexing *X,
                           int incX,
                           complexing *Y,
                           int incY) nogil
cdef complexing fused_dotc(int N,
                           complexing *X,
                           int incX,
                           complexing *Y,
                           int incY) nogil
cdef complexing fused_dotu(int N,
                           complexing *X,
                           int incX,
                           complexing *Y,
                           int incY) nogil
cdef void fused_geru(CBLAS_ORDER Order,
                     int M,
                     int N,
                     complexing alpha,
                     complexing *X,
                     int incX,
                     complexing *Y,
                     int incY,
                     complexing *A,
                     int lda) nogil
cdef double fused_nrm2(int N,
                       complexing *X,
                       int incX) nogil
