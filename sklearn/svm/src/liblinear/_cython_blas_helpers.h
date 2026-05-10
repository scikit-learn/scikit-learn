#ifndef _CYTHON_BLAS_HELPERS_H
#define _CYTHON_BLAS_HELPERS_H

typedef double (*dot_func)(int, const double*, int, const double*, int);
typedef void (*axpy_func)(int, double, const double*, int, double*, int);
typedef void (*scal_func)(int, double, const double*, int);
typedef double (*nrm2_func)(int, const double*, int);

typedef struct BlasFunctions{
    dot_func dot;
    axpy_func axpy;
    scal_func scal;
    nrm2_func nrm2;
} BlasFunctions;

#endif
