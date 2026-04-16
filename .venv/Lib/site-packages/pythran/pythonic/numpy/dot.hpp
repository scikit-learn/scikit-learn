#ifndef PYTHONIC_NUMPY_DOT_HPP
#define PYTHONIC_NUMPY_DOT_HPP

#include "pythonic/include/numpy/dot.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/multiply.hpp"
#include "pythonic/numpy/sum.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/traits.hpp"

#ifdef PYTHRAN_BLAS_NONE
#error pythran configured without BLAS but BLAS seem needed
#endif

#if defined(PYTHRAN_BLAS_SCIPY_OPENBLAS)
#define BLAS_MANGLE(s) scipy_##s##64_

/* FIXED VENDORED HEADER { */
#include "openblas_config.h"
#include <stddef.h>

extern "C" {
/* Assume C declarations for C++ */

/*Set the number of threads on runtime.*/
void scipy_openblas_set_num_threads64_(int num_threads);
void scipy_goto_set_num_threads64_(int num_threads);
int scipy_openblas_set_num_threads_local64_(int num_threads);

/*Get the number of threads on runtime.*/
int scipy_openblas_get_num_threads64_(void);

/*Get the number of physical processors (cores).*/
int scipy_openblas_get_num_procs64_(void);

/*Get the build configure on runtime.*/
char *scipy_openblas_get_config64_(void);

/*Get the CPU corename on runtime.*/
char *scipy_openblas_get_corename64_(void);

/*Set the threading backend to a custom callback.*/
typedef void (*scipy_openblas_dojob_callback64_)(int thread_num, void *jobdata, int dojob_data);
typedef void (*scipy_openblas_threads_callback64_)(int sync, scipy_openblas_dojob_callback64_ dojob,
                                                   int numjobs, size_t jobdata_elsize,
                                                   void *jobdata, int dojob_data);
void scipy_openblas_set_threads_callback_function64_(scipy_openblas_threads_callback64_ callback);

#ifdef OPENBLAS_OS_LINUX
/* Sets thread affinity for OpenBLAS threads. `thread_idx` is in [0,
 * scipy_openblas_get_num_threads64_()-1]. */
int scipy_openblas_setaffinity64_(int thread_idx, size_t cpusetsize, cpu_set_t *cpu_set);
/* Queries thread affinity for OpenBLAS threads. `thread_idx` is in [0,
 * scipy_openblas_get_num_threads64_()-1]. */
int scipy_openblas_getaffinity64_(int thread_idx, size_t cpusetsize, cpu_set_t *cpu_set);
#endif

/* Get the parallelization type which is used by OpenBLAS */
int scipy_openblas_get_parallel64_(void);
/* OpenBLAS is compiled for sequential use  */
#define OPENBLAS_SEQUENTIAL 0
/* OpenBLAS is compiled using normal threading model */
#define OPENBLAS_THREAD 1
/* OpenBLAS is compiled using OpenMP threading model */
#define OPENBLAS_OPENMP 2

/*
 * Since all of GotoBlas was written without const,
 * we disable it at build time.
 */
#ifndef OPENBLAS_CONST
#define OPENBLAS_CONST const
#endif

#define CBLAS_INDEX size_t

typedef enum CBLAS_ORDER { CblasRowMajor = 101, CblasColMajor = 102 } CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE {
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113,
  CblasConjNoTrans = 114
} CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO { CblasUpper = 121, CblasLower = 122 } CBLAS_UPLO;
typedef enum CBLAS_DIAG { CblasNonUnit = 131, CblasUnit = 132 } CBLAS_DIAG;
typedef enum CBLAS_SIDE { CblasLeft = 141, CblasRight = 142 } CBLAS_SIDE;
typedef CBLAS_ORDER CBLAS_LAYOUT;

float scipy_cblas_sdsdot64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float alpha,
                            OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx,
                            OPENBLAS_CONST float *y, OPENBLAS_CONST blasint incy);
double scipy_cblas_dsdot64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                            OPENBLAS_CONST blasint incx, OPENBLAS_CONST float *y,
                            OPENBLAS_CONST blasint incy);
float scipy_cblas_sdot64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                          OPENBLAS_CONST blasint incx, OPENBLAS_CONST float *y,
                          OPENBLAS_CONST blasint incy);
double scipy_cblas_ddot64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                           OPENBLAS_CONST blasint incx, OPENBLAS_CONST double *y,
                           OPENBLAS_CONST blasint incy);

openblas_complex_float scipy_cblas_cdotu64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                            OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                                            OPENBLAS_CONST blasint incy);
openblas_complex_float scipy_cblas_cdotc64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                            OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                                            OPENBLAS_CONST blasint incy);
openblas_complex_double scipy_cblas_zdotu64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                             OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                                             OPENBLAS_CONST blasint incy);
openblas_complex_double scipy_cblas_zdotc64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                             OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                                             OPENBLAS_CONST blasint incy);

void scipy_cblas_cdotu_sub64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                              OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                              OPENBLAS_CONST blasint incy, void *ret);
void scipy_cblas_cdotc_sub64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                              OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                              OPENBLAS_CONST blasint incy, void *ret);
void scipy_cblas_zdotu_sub64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                              OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                              OPENBLAS_CONST blasint incy, void *ret);
void scipy_cblas_zdotc_sub64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                              OPENBLAS_CONST blasint incx, OPENBLAS_CONST void *y,
                              OPENBLAS_CONST blasint incy, void *ret);

float scipy_cblas_sasum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                           OPENBLAS_CONST blasint incx);
double scipy_cblas_dasum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                            OPENBLAS_CONST blasint incx);
float scipy_cblas_scasum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                            OPENBLAS_CONST blasint incx);
double scipy_cblas_dzasum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                             OPENBLAS_CONST blasint incx);

float scipy_cblas_ssum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                          OPENBLAS_CONST blasint incx);
double scipy_cblas_dsum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                           OPENBLAS_CONST blasint incx);
float scipy_cblas_scsum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                           OPENBLAS_CONST blasint incx);
double scipy_cblas_dzsum64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                            OPENBLAS_CONST blasint incx);

float scipy_cblas_snrm264_(OPENBLAS_CONST blasint N, OPENBLAS_CONST float *X,
                           OPENBLAS_CONST blasint incX);
double scipy_cblas_dnrm264_(OPENBLAS_CONST blasint N, OPENBLAS_CONST double *X,
                            OPENBLAS_CONST blasint incX);
float scipy_cblas_scnrm264_(OPENBLAS_CONST blasint N, OPENBLAS_CONST void *X,
                            OPENBLAS_CONST blasint incX);
double scipy_cblas_dznrm264_(OPENBLAS_CONST blasint N, OPENBLAS_CONST void *X,
                             OPENBLAS_CONST blasint incX);

CBLAS_INDEX scipy_cblas_isamax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                                  OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_idamax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                                  OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_icamax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                  OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_izamax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                  OPENBLAS_CONST blasint incx);

CBLAS_INDEX scipy_cblas_isamin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                                  OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_idamin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                                  OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_icamin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                  OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_izamin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                  OPENBLAS_CONST blasint incx);

float scipy_cblas_samax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                           OPENBLAS_CONST blasint incx);
double scipy_cblas_damax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                            OPENBLAS_CONST blasint incx);
float scipy_cblas_scamax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                            OPENBLAS_CONST blasint incx);
double scipy_cblas_dzamax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                             OPENBLAS_CONST blasint incx);

float scipy_cblas_samin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                           OPENBLAS_CONST blasint incx);
double scipy_cblas_damin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                            OPENBLAS_CONST blasint incx);
float scipy_cblas_scamin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                            OPENBLAS_CONST blasint incx);
double scipy_cblas_dzamin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                             OPENBLAS_CONST blasint incx);

CBLAS_INDEX scipy_cblas_ismax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                                 OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_idmax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                                 OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_icmax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                 OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_izmax64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                 OPENBLAS_CONST blasint incx);

CBLAS_INDEX scipy_cblas_ismin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                                 OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_idmin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                                 OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_icmin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                 OPENBLAS_CONST blasint incx);
CBLAS_INDEX scipy_cblas_izmin64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                                 OPENBLAS_CONST blasint incx);

void scipy_cblas_saxpy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx, float *y,
                          OPENBLAS_CONST blasint incy);
void scipy_cblas_daxpy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx, double *y,
                          OPENBLAS_CONST blasint incy);
void scipy_cblas_caxpy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx, void *y,
                          OPENBLAS_CONST blasint incy);
void scipy_cblas_zaxpy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx, void *y,
                          OPENBLAS_CONST blasint incy);

void scipy_cblas_caxpyc64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx, void *y,
                           OPENBLAS_CONST blasint incy);
void scipy_cblas_zaxpyc64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx, void *y,
                           OPENBLAS_CONST blasint incy);

void scipy_cblas_scopy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x,
                          OPENBLAS_CONST blasint incx, float *y, OPENBLAS_CONST blasint incy);
void scipy_cblas_dcopy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x,
                          OPENBLAS_CONST blasint incx, double *y, OPENBLAS_CONST blasint incy);
void scipy_cblas_ccopy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                          OPENBLAS_CONST blasint incx, void *y, OPENBLAS_CONST blasint incy);
void scipy_cblas_zcopy64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                          OPENBLAS_CONST blasint incx, void *y, OPENBLAS_CONST blasint incy);

void scipy_cblas_sswap64_(OPENBLAS_CONST blasint n, float *x, OPENBLAS_CONST blasint incx, float *y,
                          OPENBLAS_CONST blasint incy);
void scipy_cblas_dswap64_(OPENBLAS_CONST blasint n, double *x, OPENBLAS_CONST blasint incx,
                          double *y, OPENBLAS_CONST blasint incy);
void scipy_cblas_cswap64_(OPENBLAS_CONST blasint n, void *x, OPENBLAS_CONST blasint incx, void *y,
                          OPENBLAS_CONST blasint incy);
void scipy_cblas_zswap64_(OPENBLAS_CONST blasint n, void *x, OPENBLAS_CONST blasint incx, void *y,
                          OPENBLAS_CONST blasint incy);

void scipy_cblas_srot64_(OPENBLAS_CONST blasint N, float *X, OPENBLAS_CONST blasint incX, float *Y,
                         OPENBLAS_CONST blasint incY, OPENBLAS_CONST float c,
                         OPENBLAS_CONST float s);
void scipy_cblas_drot64_(OPENBLAS_CONST blasint N, double *X, OPENBLAS_CONST blasint incX,
                         double *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST double c,
                         OPENBLAS_CONST double s);
void scipy_cblas_csrot64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                          OPENBLAS_CONST blasint incx, void *y, OPENBLAS_CONST blasint incY,
                          OPENBLAS_CONST float c, OPENBLAS_CONST float s);
void scipy_cblas_zdrot64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *x,
                          OPENBLAS_CONST blasint incx, void *y, OPENBLAS_CONST blasint incY,
                          OPENBLAS_CONST double c, OPENBLAS_CONST double s);

void scipy_cblas_srotg64_(float *a, float *b, float *c, float *s);
void scipy_cblas_drotg64_(double *a, double *b, double *c, double *s);
void scipy_cblas_crotg64_(void *a, void *b, float *c, void *s);
void scipy_cblas_zrotg64_(void *a, void *b, double *c, void *s);

void scipy_cblas_srotm64_(OPENBLAS_CONST blasint N, float *X, OPENBLAS_CONST blasint incX, float *Y,
                          OPENBLAS_CONST blasint incY, OPENBLAS_CONST float *P);
void scipy_cblas_drotm64_(OPENBLAS_CONST blasint N, double *X, OPENBLAS_CONST blasint incX,
                          double *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST double *P);

void scipy_cblas_srotmg64_(float *d1, float *d2, float *b1, OPENBLAS_CONST float b2, float *P);
void scipy_cblas_drotmg64_(double *d1, double *d2, double *b1, OPENBLAS_CONST double b2, double *P);

void scipy_cblas_sscal64_(OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, float *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_dscal64_(OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, double *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_cscal64_(OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha, void *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_zscal64_(OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha, void *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_csscal64_(OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, void *X,
                           OPENBLAS_CONST blasint incX);
void scipy_cblas_zdscal64_(OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, void *X,
                           OPENBLAS_CONST blasint incX);

void scipy_cblas_sgemv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m,
                          OPENBLAS_CONST blasint n, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST float *a, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx,
                          OPENBLAS_CONST float beta, float *y, OPENBLAS_CONST blasint incy);
void scipy_cblas_dgemv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m,
                          OPENBLAS_CONST blasint n, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST double *a, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx,
                          OPENBLAS_CONST double beta, double *y, OPENBLAS_CONST blasint incy);
void scipy_cblas_cgemv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m,
                          OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *a, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx,
                          OPENBLAS_CONST void *beta, void *y, OPENBLAS_CONST blasint incy);
void scipy_cblas_zgemv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m,
                          OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *a, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx,
                          OPENBLAS_CONST void *beta, void *y, OPENBLAS_CONST blasint incy);

void scipy_cblas_sger64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha,
                         OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX,
                         OPENBLAS_CONST float *Y, OPENBLAS_CONST blasint incY, float *A,
                         OPENBLAS_CONST blasint lda);
void scipy_cblas_dger64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha,
                         OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX,
                         OPENBLAS_CONST double *Y, OPENBLAS_CONST blasint incY, double *A,
                         OPENBLAS_CONST blasint lda);
void scipy_cblas_cgeru64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *Y, OPENBLAS_CONST blasint incY, void *A,
                          OPENBLAS_CONST blasint lda);
void scipy_cblas_cgerc64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *Y, OPENBLAS_CONST blasint incY, void *A,
                          OPENBLAS_CONST blasint lda);
void scipy_cblas_zgeru64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *Y, OPENBLAS_CONST blasint incY, void *A,
                          OPENBLAS_CONST blasint lda);
void scipy_cblas_zgerc64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *Y, OPENBLAS_CONST blasint incY, void *A,
                          OPENBLAS_CONST blasint lda);

void scipy_cblas_strsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_dtrsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_ctrsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_ztrsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *X,
                          OPENBLAS_CONST blasint incX);

void scipy_cblas_strmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_dtrmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_ctrmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *X,
                          OPENBLAS_CONST blasint incX);
void scipy_cblas_ztrmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *X,
                          OPENBLAS_CONST blasint incX);

void scipy_cblas_ssyr64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha,
                         OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, float *A,
                         OPENBLAS_CONST blasint lda);
void scipy_cblas_dsyr64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha,
                         OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, double *A,
                         OPENBLAS_CONST blasint lda);
void scipy_cblas_cher64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha,
                         OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX, void *A,
                         OPENBLAS_CONST blasint lda);
void scipy_cblas_zher64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha,
                         OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX, void *A,
                         OPENBLAS_CONST blasint lda);

void scipy_cblas_ssyr264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float alpha, OPENBLAS_CONST float *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST float *Y,
                          OPENBLAS_CONST blasint incY, float *A, OPENBLAS_CONST blasint lda);
void scipy_cblas_dsyr264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double alpha, OPENBLAS_CONST double *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST double *Y,
                          OPENBLAS_CONST blasint incY, double *A, OPENBLAS_CONST blasint lda);
void scipy_cblas_cher264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST void *Y,
                          OPENBLAS_CONST blasint incY, void *A, OPENBLAS_CONST blasint lda);
void scipy_cblas_zher264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST void *Y,
                          OPENBLAS_CONST blasint incY, void *A, OPENBLAS_CONST blasint lda);

void scipy_cblas_sgbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint KL,
                          OPENBLAS_CONST blasint KU, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST float beta, float *Y, OPENBLAS_CONST blasint incY);
void scipy_cblas_dgbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint KL,
                          OPENBLAS_CONST blasint KU, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST double beta, double *Y, OPENBLAS_CONST blasint incY);
void scipy_cblas_cgbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint KL,
                          OPENBLAS_CONST blasint KU, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *beta, void *Y, OPENBLAS_CONST blasint incY);
void scipy_cblas_zgbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint KL,
                          OPENBLAS_CONST blasint KU, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *beta, void *Y, OPENBLAS_CONST blasint incY);

void scipy_cblas_ssbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST float beta, float *Y, OPENBLAS_CONST blasint incY);
void scipy_cblas_dsbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST double beta, double *Y, OPENBLAS_CONST blasint incY);

void scipy_cblas_stbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST float *A,
                          OPENBLAS_CONST blasint lda, float *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_dtbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST double *A,
                          OPENBLAS_CONST blasint lda, double *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ctbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, void *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ztbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, void *X, OPENBLAS_CONST blasint incX);

void scipy_cblas_stbsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST float *A,
                          OPENBLAS_CONST blasint lda, float *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_dtbsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST double *A,
                          OPENBLAS_CONST blasint lda, double *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ctbsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, void *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ztbsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, void *X, OPENBLAS_CONST blasint incX);

void scipy_cblas_stpmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float *Ap, float *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_dtpmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double *Ap, double *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ctpmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *Ap, void *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ztpmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *Ap, void *X, OPENBLAS_CONST blasint incX);

void scipy_cblas_stpsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float *Ap, float *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_dtpsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double *Ap, double *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ctpsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *Ap, void *X, OPENBLAS_CONST blasint incX);
void scipy_cblas_ztpsv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *Ap, void *X, OPENBLAS_CONST blasint incX);

void scipy_cblas_ssymv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST float beta, float *Y,
                          OPENBLAS_CONST blasint incY);
void scipy_cblas_dsymv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST double beta, double *Y,
                          OPENBLAS_CONST blasint incY);
void scipy_cblas_chemv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST void *beta, void *Y,
                          OPENBLAS_CONST blasint incY);
void scipy_cblas_zhemv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST void *beta, void *Y,
                          OPENBLAS_CONST blasint incY);

void scipy_cblas_sspmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float alpha, OPENBLAS_CONST float *Ap,
                          OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST float beta, float *Y, OPENBLAS_CONST blasint incY);
void scipy_cblas_dspmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double alpha, OPENBLAS_CONST double *Ap,
                          OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST double beta, double *Y, OPENBLAS_CONST blasint incY);

void scipy_cblas_sspr64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha,
                         OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, float *Ap);
void scipy_cblas_dspr64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha,
                         OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, double *Ap);

void scipy_cblas_chpr64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha,
                         OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX, void *A);
void scipy_cblas_zhpr64_(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                         OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha,
                         OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX, void *A);

void scipy_cblas_sspr264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float alpha, OPENBLAS_CONST float *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST float *Y,
                          OPENBLAS_CONST blasint incY, float *A);
void scipy_cblas_dspr264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double alpha, OPENBLAS_CONST double *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST double *Y,
                          OPENBLAS_CONST blasint incY, double *A);
void scipy_cblas_chpr264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST void *Y,
                          OPENBLAS_CONST blasint incY, void *Ap);
void scipy_cblas_zhpr264_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *X,
                          OPENBLAS_CONST blasint incX, OPENBLAS_CONST void *Y,
                          OPENBLAS_CONST blasint incY, void *Ap);

void scipy_cblas_chbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *beta, void *Y, OPENBLAS_CONST blasint incY);
void scipy_cblas_zhbmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *beta, void *Y, OPENBLAS_CONST blasint incY);

void scipy_cblas_chpmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *Ap,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *beta, void *Y, OPENBLAS_CONST blasint incY);
void scipy_cblas_zhpmv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *Ap,
                          OPENBLAS_CONST void *X, OPENBLAS_CONST blasint incX,
                          OPENBLAS_CONST void *beta, void *Y, OPENBLAS_CONST blasint incY);

void scipy_cblas_sgemm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
                          OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta, float *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_dgemm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
                          OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta, double *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_cgemm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_cgemm3m64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                            OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                            OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                            OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
                            OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                            OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                            OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                            OPENBLAS_CONST blasint ldc);
void scipy_cblas_zgemm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_zgemm3m64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                            OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                            OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                            OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
                            OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                            OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                            OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                            OPENBLAS_CONST blasint ldc);

void scipy_cblas_sgemmt64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha,
                           OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_dgemmt64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha,
                           OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_cgemmt64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST void *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST void *beta, void *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_zgemmt64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST void *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST void *beta, void *C, OPENBLAS_CONST blasint ldc);

void scipy_cblas_ssymm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta, float *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_dsymm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta, double *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_csymm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_zsymm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                          OPENBLAS_CONST blasint ldc);

void scipy_cblas_ssyrk64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_dsyrk64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_csyrk64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *beta, void *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_zsyrk64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST void *beta, void *C, OPENBLAS_CONST blasint ldc);

void scipy_cblas_ssyr2k64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha,
                           OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_dsyr2k64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha,
                           OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_csyr2k64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST void *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST void *beta, void *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_zsyr2k64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST void *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST void *beta, void *C, OPENBLAS_CONST blasint ldc);

void scipy_cblas_strmm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *B,
                          OPENBLAS_CONST blasint ldb);
void scipy_cblas_dtrmm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *B,
                          OPENBLAS_CONST blasint ldb);
void scipy_cblas_ctrmm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *B,
                          OPENBLAS_CONST blasint ldb);
void scipy_cblas_ztrmm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *B,
                          OPENBLAS_CONST blasint ldb);

void scipy_cblas_strsm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *B,
                          OPENBLAS_CONST blasint ldb);
void scipy_cblas_dtrsm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *B,
                          OPENBLAS_CONST blasint ldb);
void scipy_cblas_ctrsm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *B,
                          OPENBLAS_CONST blasint ldb);
void scipy_cblas_ztrsm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                          OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M,
                          OPENBLAS_CONST blasint N, OPENBLAS_CONST void *alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda, void *B,
                          OPENBLAS_CONST blasint ldb);

void scipy_cblas_chemm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                          OPENBLAS_CONST blasint ldc);
void scipy_cblas_zhemm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST void *alpha, OPENBLAS_CONST void *A,
                          OPENBLAS_CONST blasint lda, OPENBLAS_CONST void *B,
                          OPENBLAS_CONST blasint ldb, OPENBLAS_CONST void *beta, void *C,
                          OPENBLAS_CONST blasint ldc);

void scipy_cblas_cherk64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST float beta, void *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_zherk64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                          OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                          OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                          OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha,
                          OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                          OPENBLAS_CONST double beta, void *C, OPENBLAS_CONST blasint ldc);

void scipy_cblas_cher2k64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST void *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST float beta, void *C, OPENBLAS_CONST blasint ldc);
void scipy_cblas_zher2k64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                           OPENBLAS_CONST blasint K, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *A, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST void *B, OPENBLAS_CONST blasint ldb,
                           OPENBLAS_CONST double beta, void *C, OPENBLAS_CONST blasint ldc);

void scipy_cblas_xerbla64_(blasint p, OPENBLAS_CONST char *rout, OPENBLAS_CONST char *form, ...);

/*** BLAS extensions ***/

void scipy_cblas_saxpby64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float alpha,
                           OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx,
                           OPENBLAS_CONST float beta, float *y, OPENBLAS_CONST blasint incy);

void scipy_cblas_daxpby64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double alpha,
                           OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx,
                           OPENBLAS_CONST double beta, double *y, OPENBLAS_CONST blasint incy);

void scipy_cblas_caxpby64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx,
                           OPENBLAS_CONST void *beta, void *y, OPENBLAS_CONST blasint incy);

void scipy_cblas_zaxpby64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST void *alpha,
                           OPENBLAS_CONST void *x, OPENBLAS_CONST blasint incx,
                           OPENBLAS_CONST void *beta, void *y, OPENBLAS_CONST blasint incy);

void scipy_cblas_somatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST float calpha, OPENBLAS_CONST float *a,
                              OPENBLAS_CONST blasint clda, float *b, OPENBLAS_CONST blasint cldb);
void scipy_cblas_domatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST double calpha, OPENBLAS_CONST double *a,
                              OPENBLAS_CONST blasint clda, double *b, OPENBLAS_CONST blasint cldb);
void scipy_cblas_comatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST float *calpha, OPENBLAS_CONST float *a,
                              OPENBLAS_CONST blasint clda, float *b, OPENBLAS_CONST blasint cldb);
void scipy_cblas_zomatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST double *calpha, OPENBLAS_CONST double *a,
                              OPENBLAS_CONST blasint clda, double *b, OPENBLAS_CONST blasint cldb);

void scipy_cblas_simatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST float calpha, float *a, OPENBLAS_CONST blasint clda,
                              OPENBLAS_CONST blasint cldb);
void scipy_cblas_dimatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST double calpha, double *a, OPENBLAS_CONST blasint clda,
                              OPENBLAS_CONST blasint cldb);
void scipy_cblas_cimatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST float *calpha, float *a, OPENBLAS_CONST blasint clda,
                              OPENBLAS_CONST blasint cldb);
void scipy_cblas_zimatcopy64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
                              OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
                              OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
                              OPENBLAS_CONST double *calpha, double *a, OPENBLAS_CONST blasint clda,
                              OPENBLAS_CONST blasint cldb);

void scipy_cblas_sgeadd64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST blasint crows,
                           OPENBLAS_CONST blasint ccols, OPENBLAS_CONST float calpha, float *a,
                           OPENBLAS_CONST blasint clda, OPENBLAS_CONST float cbeta, float *c,
                           OPENBLAS_CONST blasint cldc);
void scipy_cblas_dgeadd64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST blasint crows,
                           OPENBLAS_CONST blasint ccols, OPENBLAS_CONST double calpha, double *a,
                           OPENBLAS_CONST blasint clda, OPENBLAS_CONST double cbeta, double *c,
                           OPENBLAS_CONST blasint cldc);
void scipy_cblas_cgeadd64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST blasint crows,
                           OPENBLAS_CONST blasint ccols, OPENBLAS_CONST float *calpha, float *a,
                           OPENBLAS_CONST blasint clda, OPENBLAS_CONST float *cbeta, float *c,
                           OPENBLAS_CONST blasint cldc);
void scipy_cblas_zgeadd64_(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST blasint crows,
                           OPENBLAS_CONST blasint ccols, OPENBLAS_CONST double *calpha, double *a,
                           OPENBLAS_CONST blasint clda, OPENBLAS_CONST double *cbeta, double *c,
                           OPENBLAS_CONST blasint cldc);

/*** BFLOAT16 and INT8 extensions ***/
/* convert float array to BFLOAT16 array by rounding */
void scipy_cblas_sbstobf1664_(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *in,
                              OPENBLAS_CONST blasint incin, bfloat16 *out,
                              OPENBLAS_CONST blasint incout);
/* convert double array to BFLOAT16 array by rounding */
void scipy_cblas_sbdtobf1664_(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *in,
                              OPENBLAS_CONST blasint incin, bfloat16 *out,
                              OPENBLAS_CONST blasint incout);
/* convert BFLOAT16 array to float array */
void scipy_cblas_sbf16tos64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST bfloat16 *in,
                             OPENBLAS_CONST blasint incin, float *out,
                             OPENBLAS_CONST blasint incout);
/* convert BFLOAT16 array to double array */
void scipy_cblas_dbf16tod64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST bfloat16 *in,
                             OPENBLAS_CONST blasint incin, double *out,
                             OPENBLAS_CONST blasint incout);
/* dot production of BFLOAT16 input arrays, and output as float */
float scipy_cblas_sbdot64_(OPENBLAS_CONST blasint n, OPENBLAS_CONST bfloat16 *x,
                           OPENBLAS_CONST blasint incx, OPENBLAS_CONST bfloat16 *y,
                           OPENBLAS_CONST blasint incy);
void scipy_cblas_sbgemv64_(OPENBLAS_CONST enum CBLAS_ORDER order,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m,
                           OPENBLAS_CONST blasint n, OPENBLAS_CONST float alpha,
                           OPENBLAS_CONST bfloat16 *a, OPENBLAS_CONST blasint lda,
                           OPENBLAS_CONST bfloat16 *x, OPENBLAS_CONST blasint incx,
                           OPENBLAS_CONST float beta, float *y, OPENBLAS_CONST blasint incy);

void scipy_cblas_sbgemm64_(OPENBLAS_CONST enum CBLAS_ORDER Order,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                           OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M,
                           OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
                           OPENBLAS_CONST float alpha, OPENBLAS_CONST bfloat16 *A,
                           OPENBLAS_CONST blasint lda, OPENBLAS_CONST bfloat16 *B,
                           OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta, float *C,
                           OPENBLAS_CONST blasint ldc);
}
/* } VENDORED HEADER */

#else
#define BLAS_MANGLE(s) s
#if defined(PYTHRAN_BLAS_ATLAS) || defined(PYTHRAN_BLAS_SATLAS)
extern "C" {
#endif
#include <cblas.h>
#if defined(PYTHRAN_BLAS_ATLAS) || defined(PYTHRAN_BLAS_SATLAS)
}
#endif
#endif

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  std::enable_if_t<types::is_dtype<E>::value && types::is_dtype<F>::value,
                   decltype(std::declval<E>() * std::declval<F>())>
  dot(E const &e, F const &f)
  {
    return e * f;
  }

  template <class E>
  struct blas_buffer_t {
    typename E::dtype const *operator()(E const &e) const
    {
      return e.buffer;
    }
  };
  template <class T>
  struct blas_buffer_t<types::list<T>> {
    T const *operator()(types::list<T> const &e) const
    {
      return &e.fast(0);
    }
  };
  template <class T, size_t N>
  struct blas_buffer_t<types::array_tuple<T, N>> {
    T const *operator()(types::array_tuple<T, N> const &e) const
    {
      return e.data();
    }
  };

  template <class E>
  auto blas_buffer(E const &e) -> decltype(blas_buffer_t<E>{}(e))
  {
    return blas_buffer_t<E>{}(e);
  }

  template <class E, class... S>
  typename E::dtype const *blas_buffer(types::numpy_gexpr<E, S...> const &e)
  {
    return e.data();
  }

  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value &&
                       types::is_numexpr_arg<F>::value   // Arguments are array_like
                       && E::value == 1 && F::value == 1 // It is a two vectors.
                       && (!is_blas_view<E>::value || !is_blas_view<F>::value ||
                           !std::is_same<typename E::dtype, typename F::dtype>::value),
                   typename __combined<typename E::dtype, typename F::dtype>::type>
  dot(E const &e, F const &f)
  {
    return sum(functor::multiply{}(e, f));
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, float>::value &&
                       std::is_same<typename F::dtype, float>::value && is_blas_array<E>::value &&
                       is_blas_array<F>::value,
                   float>
  dot(E const &e, F const &f)
  {
    return BLAS_MANGLE(cblas_sdot)(e.size(), blas_buffer(e), 1, blas_buffer(f), 1);
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, double>::value &&
                       std::is_same<typename F::dtype, double>::value && is_blas_array<E>::value &&
                       is_blas_array<F>::value,
                   double>
  dot(E const &e, F const &f)
  {
    return BLAS_MANGLE(cblas_ddot)(e.size(), blas_buffer(e), 1, blas_buffer(f), 1);
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<float>>::value &&
                       std::is_same<typename F::dtype, std::complex<float>>::value &&
                       is_blas_array<E>::value && is_blas_array<F>::value,
                   std::complex<float>>
  dot(E const &e, F const &f)
  {
    std::complex<float> out;
    BLAS_MANGLE(cblas_cdotu_sub)
    (e.size(), blas_buffer(e), 1, blas_buffer(f), 1, &out);
    return out;
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<double>>::value &&
                       std::is_same<typename F::dtype, std::complex<double>>::value &&
                       is_blas_array<E>::value && is_blas_array<F>::value,
                   std::complex<double>>
  dot(E const &e, F const &f)
  {
    std::complex<double> out;
    BLAS_MANGLE(cblas_zdotu_sub)
    (e.size(), blas_buffer(e), 1, blas_buffer(f), 1, &out);
    return out;
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, float>::value &&
                       std::is_same<typename F::dtype, float>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   float>
  dot(E const &e, F const &f)
  {
    if (e.template strides<0>() >= 1 && f.template strides<0>() >= 1) {
      return BLAS_MANGLE(cblas_sdot)(e.size(), blas_buffer(e), e.template strides<0>(),
                                     blas_buffer(f), f.template strides<0>());
    } else {
      return dot(asarray(e), asarray(f));
    }
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, double>::value &&
                       std::is_same<typename F::dtype, double>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   double>
  dot(E const &e, F const &f)
  {
    if (e.template strides<0>() >= 1 && f.template strides<0>() >= 1) {
      return BLAS_MANGLE(cblas_ddot)(e.size(), blas_buffer(e), e.template strides<0>(),
                                     blas_buffer(f), f.template strides<0>());
    } else {
      return dot(asarray(e), asarray(f));
    }
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<float>>::value &&
                       std::is_same<typename F::dtype, std::complex<float>>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   std::complex<float>>
  dot(E const &e, F const &f)
  {
    if (e.template strides<0>() >= 1 && f.template strides<0>() >= 1) {
      std::complex<float> out;
      BLAS_MANGLE(cblas_cdotu_sub)
      (e.size(), blas_buffer(e), e.template strides<0>(), blas_buffer(f), f.template strides<0>(),
       &out);
      return out;
    } else {
      return dot(asarray(e), asarray(f));
    }
  }

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<double>>::value &&
                       std::is_same<typename F::dtype, std::complex<double>>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   std::complex<double>>
  dot(E const &e, F const &f)
  {
    if (e.template strides<0>() >= 1 && f.template strides<0>() >= 1) {
      std::complex<double> out;
      BLAS_MANGLE(cblas_zdotu_sub)
      (e.size(), blas_buffer(e), e.template strides<0>(), blas_buffer(f), f.template strides<0>(),
       &out);
      return out;
    } else {
      return dot(asarray(e), asarray(f));
    }
  }

  /// Matrice / Vector multiplication

#define MV_DEF(T, L)                                                                               \
  inline void mv(int m, int n, T *A, T *B, T *C)                                                   \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemv)                                                                   \
    (CblasRowMajor, CblasNoTrans, n, m, 1, A, m, B, 1, 0, C, 1);                                   \
  }

  MV_DEF(double, d)
  MV_DEF(float, s)

#undef MV_DEF

#define TV_DEF(T, L)                                                                               \
  inline void tv(int m, int n, T *A, T *B, T *C)                                                   \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemv)                                                                   \
    (CblasRowMajor, CblasTrans, m, n, 1, A, n, B, 1, 0, C, 1);                                     \
  }

  TV_DEF(double, d)
  TV_DEF(float, s)

#undef TV_DEF

#define MV_DEF(T, K, L)                                                                            \
  inline void mv(int m, int n, T *A, T *B, T *C)                                                   \
  {                                                                                                \
    T alpha = 1, beta = 0;                                                                         \
    BLAS_MANGLE(cblas_##L##gemv)                                                                   \
    (CblasRowMajor, CblasNoTrans, n, m, (K *)&alpha, (K *)A, m, (K *)B, 1, (K *)&beta, (K *)C, 1); \
  }
  MV_DEF(std::complex<float>, float, c)
  MV_DEF(std::complex<double>, double, z)
#undef MV_DEF

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 1,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::ndarray<E, pS0> const &f, types::ndarray<E, pS1> const &e)
  {
    types::ndarray<E, types::pshape<long>> out(types::pshape<long>{f.template shape<0>()},
                                               builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    mv(m, n, f.buffer, e.buffer, out.buffer);
    return out;
  }

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 1,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &f, types::ndarray<E, pS1> const &e)
  {
    types::ndarray<E, types::pshape<long>> out(types::pshape<long>{f.template shape<0>()},
                                               builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    tv(m, n, f.arg.buffer, e.buffer, out.buffer);
    return out;
  }

// The trick is to not transpose the matrix so that MV become VM
#define VM_DEF(T, L)                                                                               \
  inline void vm(int m, int n, T *A, T *B, T *C)                                                   \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemv)                                                                   \
    (CblasRowMajor, CblasTrans, n, m, 1, A, m, B, 1, 0, C, 1);                                     \
  }

  VM_DEF(double, d)
  VM_DEF(float, s)

#undef VM_DEF
#define VT_DEF(T, L)                                                                               \
  inline void vt(int m, int n, T *A, T *B, T *C)                                                   \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemv)                                                                   \
    (CblasRowMajor, CblasNoTrans, m, n, 1, A, n, B, 1, 0, C, 1);                                   \
  }

  VT_DEF(double, d)
  VT_DEF(float, s)

#undef VM_DEF
#define VM_DEF(T, K, L)                                                                            \
  inline void vm(int m, int n, T *A, T *B, T *C)                                                   \
  {                                                                                                \
    T alpha = 1, beta = 0;                                                                         \
    BLAS_MANGLE(cblas_##L##gemv)                                                                   \
    (CblasRowMajor, CblasTrans, n, m, (K *)&alpha, (K *)A, m, (K *)B, 1, (K *)&beta, (K *)C, 1);   \
  }
  VM_DEF(std::complex<float>, float, c)
  VM_DEF(std::complex<double>, double, z)
#undef VM_DEF

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 1 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::ndarray<E, pS0> const &e, types::ndarray<E, pS1> const &f)
  {
    types::ndarray<E, types::pshape<long>> out(types::pshape<long>{f.template shape<1>()},
                                               builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    vm(m, n, f.buffer, e.buffer, out.buffer);
    return out;
  }

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 1 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::ndarray<E, pS0> const &e, types::numpy_texpr<types::ndarray<E, pS1>> const &f)
  {
    types::ndarray<E, types::pshape<long>> out(types::pshape<long>{f.template shape<1>()},
                                               builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    vt(m, n, f.arg.buffer, e.buffer, out.buffer);
    return out;
  }

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value &&
                       types::is_numexpr_arg<F>::value // It is an array_like
                       && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
                           !std::is_same<typename E::dtype, typename F::dtype>::value) &&
                       is_blas_type<typename E::dtype>::value &&
                       is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                              // blas
                       && E::value == 2 && F::value == 1,     // And it is matrix / vect
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f)
  {
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   typename E::shape_t>
        e_ = e;
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   typename F::shape_t>
        f_ = f;
    return dot(e_, f_);
  }

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value &&
                       types::is_numexpr_arg<F>::value // It is an array_like
                       && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
                           !std::is_same<typename E::dtype, typename F::dtype>::value) &&
                       is_blas_type<typename E::dtype>::value &&
                       is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                              // blas
                       && E::value == 1 && F::value == 2,     // And it is vect / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f)
  {
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   typename E::shape_t>
        e_ = e;
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   typename F::shape_t>
        f_ = f;
    return dot(e_, f_);
  }

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix vector multiplication.
  template <class E, class F>
  std::enable_if_t<(!is_blas_type<typename E::dtype>::value ||
                    !is_blas_type<typename F::dtype>::value) &&
                       E::value == 1 && F::value == 2, // And it is vect / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f)
  {
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   types::pshape<long>>
        out(types::pshape<long>{f.template shape<1>()}, 0);
    for (long i = 0; i < out.template shape<0>(); i++)
      for (long j = 0; j < f.template shape<0>(); j++)
        out[i] += e[j] * f[types::array_tuple<long, 2>{{j, i}}];
    return out;
  }

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix vector multiplication.
  template <class E, class F>
  std::enable_if_t<(!is_blas_type<typename E::dtype>::value ||
                    !is_blas_type<typename F::dtype>::value) &&
                       E::value == 2 && F::value == 1, // And it is vect / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f)
  {
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   types::pshape<long>>
        out(types::pshape<long>{e.template shape<0>()}, 0);
    for (long i = 0; i < out.template shape<0>(); i++)
      for (long j = 0; j < f.template shape<0>(); j++)
        out[i] += e[types::array_tuple<long, 2>{{i, j}}] * f[j];
    return out;
  }

  /// Matrix / Matrix multiplication

#define MM_DEF(T, L)                                                                               \
  inline void mm(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, k, B, n, 0, C, n);                  \
  }
  MM_DEF(double, d)
  MM_DEF(float, s)
#undef MM_DEF
#define MM_DEF(T, K, L)                                                                            \
  inline void mm(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    T alpha = 1, beta = 0;                                                                         \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, (K *)&alpha, (K *)A, k, (K *)B, n,        \
     (K *)&beta, (K *)C, n);                                                                       \
  }
  MM_DEF(std::complex<float>, float, c)
  MM_DEF(std::complex<double>, double, z)
#undef MM_DEF

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::ndarray<E, pS0> const &a, types::ndarray<E, pS1> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(), k = b.template shape<0>();

    types::ndarray<E, types::array_tuple<long, 2>> out(types::array_tuple<long, 2>{{m, n}},
                                                       builtins::None);
    mm(m, n, k, a.buffer, b.buffer, out.buffer);
    return out;
  }

  template <class E, class pS0, class pS1, class pS2>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2 && std::tuple_size<pS2>::value == 2,
                   types::ndarray<E, pS2>> &
  dot(types::ndarray<E, pS0> const &a, types::ndarray<E, pS1> const &b, types::ndarray<E, pS2> &c)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(), k = b.template shape<0>();

    mm(m, n, k, a.buffer, b.buffer, c.buffer);
    return c;
  }

#define TM_DEF(T, L)                                                                               \
  inline void tm(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1, A, m, B, n, 0, C, n);                    \
  }
  TM_DEF(double, d)
  TM_DEF(float, s)
#undef TM_DEF
#define TM_DEF(T, K, L)                                                                            \
  inline void tm(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    T alpha = 1, beta = 0;                                                                         \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, (K *)&alpha, (K *)A, m, (K *)B, n,          \
     (K *)&beta, (K *)C, n);                                                                       \
  }
  TM_DEF(std::complex<float>, float, c)
  TM_DEF(std::complex<double>, double, z)
#undef TM_DEF

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &a, types::ndarray<E, pS1> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(), k = b.template shape<0>();

    types::ndarray<E, types::array_tuple<long, 2>> out(types::array_tuple<long, 2>{{m, n}},
                                                       builtins::None);
    tm(m, n, k, a.arg.buffer, b.buffer, out.buffer);
    return out;
  }

#define MT_DEF(T, L)                                                                               \
  inline void mt(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1, A, k, B, k, 0, C, n);                    \
  }
  MT_DEF(double, d)
  MT_DEF(float, s)
#undef MT_DEF
#define MT_DEF(T, K, L)                                                                            \
  inline void mt(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    T alpha = 1, beta = 0;                                                                         \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, (K *)&alpha, (K *)A, k, (K *)B, k,          \
     (K *)&beta, (K *)C, n);                                                                       \
  }
  MT_DEF(std::complex<float>, float, c)
  MT_DEF(std::complex<double>, double, z)
#undef MT_DEF

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::ndarray<E, pS0> const &a, types::numpy_texpr<types::ndarray<E, pS1>> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(), k = b.template shape<0>();

    types::ndarray<E, types::array_tuple<long, 2>> out(types::array_tuple<long, 2>{{m, n}},
                                                       builtins::None);
    mt(m, n, k, a.buffer, b.arg.buffer, out.buffer);
    return out;
  }

#define TT_DEF(T, L)                                                                               \
  inline void tt(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasTrans, CblasTrans, m, n, k, 1, A, m, B, k, 0, C, n);                      \
  }
  TT_DEF(double, d)
  TT_DEF(float, s)
#undef TT_DEF
#define TT_DEF(T, K, L)                                                                            \
  inline void tt(int m, int n, int k, T *A, T *B, T *C)                                            \
  {                                                                                                \
    T alpha = 1, beta = 0;                                                                         \
    BLAS_MANGLE(cblas_##L##gemm)                                                                   \
    (CblasRowMajor, CblasTrans, CblasTrans, m, n, k, (K *)&alpha, (K *)A, m, (K *)B, k,            \
     (K *)&beta, (K *)C, n);                                                                       \
  }
  TT_DEF(std::complex<float>, float, c)
  TT_DEF(std::complex<double>, double, z)
#undef TT_DEF

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &a,
      types::numpy_texpr<types::ndarray<E, pS1>> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(), k = b.template shape<0>();

    types::ndarray<E, types::array_tuple<long, 2>> out(types::array_tuple<long, 2>{{m, n}},
                                                       builtins::None);
    tt(m, n, k, a.arg.buffer, b.arg.buffer, out.buffer);
    return out;
  }

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value &&
                       types::is_numexpr_arg<F>::value // It is an array_like
                       && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
                           !std::is_same<typename E::dtype, typename F::dtype>::value) &&
                       is_blas_type<typename E::dtype>::value &&
                       is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                              // blas
                       && E::value == 2 && F::value == 2,     // And both are matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, 2>>>
  dot(E const &e, F const &f)
  {
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   typename E::shape_t>
        e_ = e;
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   typename F::shape_t>
        f_ = f;
    return dot(e_, f_);
  }

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix multiplication.
  template <class E, class F>
  std::enable_if_t<(!is_blas_type<typename E::dtype>::value ||
                    !is_blas_type<typename F::dtype>::value) &&
                       E::value == 2 && F::value == 2, // And it is matrix / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, 2>>>
  dot(E const &e, F const &f)
  {
    types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                   types::array_tuple<long, 2>>
        out(types::array_tuple<long, 2>{{e.template shape<0>(), f.template shape<1>()}}, 0);
    for (long i = 0; i < out.template shape<0>(); i++)
      for (long j = 0; j < out.template shape<1>(); j++)
        for (long k = 0; k < e.template shape<1>(); k++)
          out[types::array_tuple<long, 2>{{i, j}}] +=
              e[types::array_tuple<long, 2>{{i, k}}] * f[types::array_tuple<long, 2>{{k, j}}];
    return out;
  }

  template <class E, class F>
  std::enable_if_t<(E::value >= 3 && F::value == 1), // And it is matrix / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, E::value - 1>>>
  dot(E const &e, F const &f)
  {
    auto out = dot(e.reshape(types::array_tuple<long, 2>{{sutils::prod_head(e), f.size()}}), f);
    types::array_tuple<long, E::value - 1> out_shape;
    auto tmp = sutils::getshape(e);
    std::copy(tmp.begin(), tmp.end() - 1, out_shape.begin());
    return out.reshape(out_shape);
  }

  template <class E, class F>
  std::enable_if_t<(E::value >= 3 && F::value >= 2),
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, E::value - 1>>>
  dot(E const &e, F const &f)
  {
    static_assert(E::value == 0, "not implemented yet");
  }
} // namespace numpy
PYTHONIC_NS_END

#undef BLAS_MANGLE

#endif
