#ifndef _SVM_CYTHON_BLAS_HELPERS_H
#define _SVM_CYTHON_BLAS_HELPERS_H

enum BLAS_Order {RowMajor=101, ColMajor=102};
enum BLAS_Trans {NoTrans=110, Trans=116};

typedef double (*dot_func)(int, double*, int, double*, int);
typedef void (*dscal_func)(int, double, double*, int);
typedef void (*dgemv_func)(enum BLAS_Order, enum BLAS_Trans, int, int, double,
                double *, int, double *, int,
                double, double *, int);

typedef struct BlasFunctions{
    dot_func dot;
    dscal_func dscal;
    dgemv_func dgemv;
} BlasFunctions;

#endif