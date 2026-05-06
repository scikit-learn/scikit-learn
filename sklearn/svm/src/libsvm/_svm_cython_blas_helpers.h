#ifndef _SVM_CYTHON_BLAS_HELPERS_H
#define _SVM_CYTHON_BLAS_HELPERS_H

typedef double (*dot_func)(int, const double*, int, const double*, int);
typedef struct BlasFunctions{
    dot_func dot;
} BlasFunctions;

#endif
