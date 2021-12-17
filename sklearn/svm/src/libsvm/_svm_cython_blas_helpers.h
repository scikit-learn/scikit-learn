#ifndef _SVM_CYTHON_BLAS_HELPERS_H
#define _SVM_CYTHON_BLAS_HELPERS_H

typedef double (*dot_func)(int, double*, int, double*, int);
typedef struct BlasFunctions{
    dot_func dot;
} BlasFunctions;

#endif
