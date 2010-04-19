/* blasp.h  --  C prototypes for BLAS                         Ver 1.0 */
/* Jesse Bennett                                       March 23, 2000 */

/* Functions  listed in alphabetical order */

#ifdef F2C_COMPAT

void cdotc_(fcomplex *dotval, int *n, fcomplex *cx, int *incx,
            fcomplex *cy, int *incy);

void cdotu_(fcomplex *dotval, int *n, fcomplex *cx, int *incx,
            fcomplex *cy, int *incy);

double sasum_(int *n, float *sx, int *incx);

double scasum_(int *n, fcomplex *cx, int *incx);

double scnrm2_(int *n, fcomplex *x, int *incx);

double sdot_(int *n, float *sx, int *incx, float *sy, int *incy);

double snrm2_(int *n, float *x, int *incx);

void zdotc_(dcomplex *dotval, int *n, dcomplex *cx, int *incx,
            dcomplex *cy, int *incy);

void zdotu_(dcomplex *dotval, int *n, dcomplex *cx, int *incx,
            dcomplex *cy, int *incy);

#else

fcomplex cdotc_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

fcomplex cdotu_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

float sasum_(int *n, float *sx, int *incx);

float scasum_(int *n, fcomplex *cx, int *incx);

float scnrm2_(int *n, fcomplex *x, int *incx);

float sdot_(int *n, float *sx, int *incx, float *sy, int *incy);

float snrm2_(int *n, float *x, int *incx);

dcomplex zdotc_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

dcomplex zdotu_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

#endif

/* Remaining functions listed in alphabetical order */

int caxpy_(int *n, fcomplex *ca, fcomplex *cx, int *incx, fcomplex *cy,
           int *incy);

int ccopy_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

int cgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           fcomplex *alpha, fcomplex *a, int *lda, fcomplex *x, int *incx,
           fcomplex *beta, fcomplex *y, int *incy);

int cgemm_(char *transa, char *transb, int *m, int *n, int *k,
           fcomplex *alpha, fcomplex *a, int *lda, fcomplex *b, int *ldb,
           fcomplex *beta, fcomplex *c, int *ldc);

int cgemv_(char *trans, int *m, int *n, fcomplex *alpha, fcomplex *a,
           int *lda, fcomplex *x, int *incx, fcomplex *beta, fcomplex *y,
           int *incy);

int cgerc_(int *m, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *a, int *lda);

int cgeru_(int *m, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *a, int *lda);

int chbmv_(char *uplo, int *n, int *k, fcomplex *alpha, fcomplex *a,
           int *lda, fcomplex *x, int *incx, fcomplex *beta, fcomplex *y,
           int *incy);

int chemm_(char *side, char *uplo, int *m, int *n, fcomplex *alpha,
           fcomplex *a, int *lda, fcomplex *b, int *ldb, fcomplex *beta,
           fcomplex *c, int *ldc);

int chemv_(char *uplo, int *n, fcomplex *alpha, fcomplex *a, int *lda,
           fcomplex *x, int *incx, fcomplex *beta, fcomplex *y, int *incy);

int cher_(char *uplo, int *n, float *alpha, fcomplex *x, int *incx,
          fcomplex *a, int *lda);

int cher2_(char *uplo, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *a, int *lda);

int cher2k_(char *uplo, char *trans, int *n, int *k, fcomplex *alpha,
            fcomplex *a, int *lda, fcomplex *b, int *ldb, float *beta,
            fcomplex *c, int *ldc);

int cherk_(char *uplo, char *trans, int *n, int *k, float *alpha,
           fcomplex *a, int *lda, float *beta, fcomplex *c, int *ldc);

int chpmv_(char *uplo, int *n, fcomplex *alpha, fcomplex *ap, fcomplex *x,
           int *incx, fcomplex *beta, fcomplex *y, int *incy);

int chpr_(char *uplo, int *n, float *alpha, fcomplex *x, int *incx,
          fcomplex *ap);

int chpr2_(char *uplo, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *ap);

int crotg_(fcomplex *ca, fcomplex *cb, float *c, fcomplex *s);

int cscal_(int *n, fcomplex *ca, fcomplex *cx, int *incx);

int csscal_(int *n, float *sa, fcomplex *cx, int *incx);

int cswap_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

int csymm_(char *side, char *uplo, int *m, int *n, fcomplex *alpha,
           fcomplex *a, int *lda, fcomplex *b, int *ldb, fcomplex *beta,
           fcomplex *c, int *ldc);

int csyr2k_(char *uplo, char *trans, int *n, int *k, fcomplex *alpha,
            fcomplex *a, int *lda, fcomplex *b, int *ldb, fcomplex *beta,
            fcomplex *c, int *ldc);

int csyrk_(char *uplo, char *trans, int *n, int *k, fcomplex *alpha,
           fcomplex *a, int *lda, fcomplex *beta, fcomplex *c, int *ldc);

int ctbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           fcomplex *a, int *lda, fcomplex *x, int *incx);

int ctbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           fcomplex *a, int *lda, fcomplex *x, int *incx);

int ctpmv_(char *uplo, char *trans, char *diag, int *n, fcomplex *ap,
           fcomplex *x, int *incx);

int ctpsv_(char *uplo, char *trans, char *diag, int *n, fcomplex *ap,
           fcomplex *x, int *incx);

int ctrmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, fcomplex *alpha, fcomplex *a, int *lda, fcomplex *b,
           int *ldb);

int ctrmv_(char *uplo, char *trans, char *diag, int *n, fcomplex *a,
           int *lda, fcomplex *x, int *incx);

int ctrsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, fcomplex *alpha, fcomplex *a, int *lda, fcomplex *b,
           int *ldb);

int ctrsv_(char *uplo, char *trans, char *diag, int *n, fcomplex *a,
           int *lda, fcomplex *x, int *incx);

int daxpy_(int *n, double *sa, double *sx, int *incx, double *sy,
           int *incy);

int dcopy_(int *n, double *sx, int *incx, double *sy, int *incy);

int dgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           double *alpha, double *a, int *lda, double *x, int *incx,
           double *beta, double *y, int *incy);

int dgemm_(char *transa, char *transb, int *m, int *n, int *k,
           double *alpha, double *a, int *lda, double *b, int *ldb,
           double *beta, double *c, int *ldc);

int dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
           int *lda, double *x, int *incx, double *beta, double *y, 
           int *incy);

int dger_(int *m, int *n, double *alpha, double *x, int *incx,
          double *y, int *incy, double *a, int *lda);

int drot_(int *n, double *sx, int *incx, double *sy, int *incy,
          double *c, double *s);

int drotg_(double *sa, double *sb, double *c, double *s);

int dsbmv_(char *uplo, int *n, int *k, double *alpha, double *a,
           int *lda, double *x, int *incx, double *beta, double *y, 
           int *incy);

int dscal_(int *n, double *sa, double *sx, int *incx);

int dspmv_(char *uplo, int *n, double *alpha, double *ap, double *x,
           int *incx, double *beta, double *y, int *incy);

int dspr_(char *uplo, int *n, double *alpha, double *x, int *incx,
          double *ap);

int dspr2_(char *uplo, int *n, double *alpha, double *x, int *incx,
           double *y, int *incy, double *ap);

int dswap_(int *n, double *sx, int *incx, double *sy, int *incy);

int dsymm_(char *side, char *uplo, int *m, int *n, double *alpha,
           double *a, int *lda, double *b, int *ldb, double *beta,
           double *c, int *ldc);

int dsymv_(char *uplo, int *n, double *alpha, double *a, int *lda,
           double *x, int *incx, double *beta, double *y, int *incy);

int dsyr_(char *uplo, int *n, double *alpha, double *x, int *incx,
          double *a, int *lda);

int dsyr2_(char *uplo, int *n, double *alpha, double *x, int *incx,
           double *y, int *incy, double *a, int *lda);

int dsyr2k_(char *uplo, char *trans, int *n, int *k, double *alpha,
            double *a, int *lda, double *b, int *ldb, double *beta,
            double *c, int *ldc);

int dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha,
           double *a, int *lda, double *beta, double *c, int *ldc);

int dtbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           double *a, int *lda, double *x, int *incx);

int dtbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           double *a, int *lda, double *x, int *incx);

int dtpmv_(char *uplo, char *trans, char *diag, int *n, double *ap,
           double *x, int *incx);

int dtpsv_(char *uplo, char *trans, char *diag, int *n, double *ap,
           double *x, int *incx);

int dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, double *alpha, double *a, int *lda, double *b, 
           int *ldb);

int dtrmv_(char *uplo, char *trans, char *diag, int *n, double *a,
           int *lda, double *x, int *incx);

int dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, double *alpha, double *a, int *lda, double *b, 
           int *ldb);

int dtrsv_(char *uplo, char *trans, char *diag, int *n, double *a,
           int *lda, double *x, int *incx);


int saxpy_(int *n, float *sa, float *sx, int *incx, float *sy, int *incy);

int scopy_(int *n, float *sx, int *incx, float *sy, int *incy);

int sgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           float *alpha, float *a, int *lda, float *x, int *incx,
           float *beta, float *y, int *incy);

int sgemm_(char *transa, char *transb, int *m, int *n, int *k,
           float *alpha, float *a, int *lda, float *b, int *ldb,
           float *beta, float *c, int *ldc);

int sgemv_(char *trans, int *m, int *n, float *alpha, float *a,
           int *lda, float *x, int *incx, float *beta, float *y, 
           int *incy);

int sger_(int *m, int *n, float *alpha, float *x, int *incx,
          float *y, int *incy, float *a, int *lda);

int srot_(int *n, float *sx, int *incx, float *sy, int *incy,
          float *c, float *s);

int srotg_(float *sa, float *sb, float *c, float *s);

int ssbmv_(char *uplo, int *n, int *k, float *alpha, float *a,
           int *lda, float *x, int *incx, float *beta, float *y, 
           int *incy);

int sscal_(int *n, float *sa, float *sx, int *incx);

int sspmv_(char *uplo, int *n, float *alpha, float *ap, float *x,
           int *incx, float *beta, float *y, int *incy);

int sspr_(char *uplo, int *n, float *alpha, float *x, int *incx,
          float *ap);

int sspr2_(char *uplo, int *n, float *alpha, float *x, int *incx,
           float *y, int *incy, float *ap);

int sswap_(int *n, float *sx, int *incx, float *sy, int *incy);

int ssymm_(char *side, char *uplo, int *m, int *n, float *alpha,
           float *a, int *lda, float *b, int *ldb, float *beta,
           float *c, int *ldc);

int ssymv_(char *uplo, int *n, float *alpha, float *a, int *lda,
           float *x, int *incx, float *beta, float *y, int *incy);

int ssyr_(char *uplo, int *n, float *alpha, float *x, int *incx,
          float *a, int *lda);

int ssyr2_(char *uplo, int *n, float *alpha, float *x, int *incx,
           float *y, int *incy, float *a, int *lda);

int ssyr2k_(char *uplo, char *trans, int *n, int *k, float *alpha,
            float *a, int *lda, float *b, int *ldb, float *beta,
            float *c, int *ldc);

int ssyrk_(char *uplo, char *trans, int *n, int *k, float *alpha,
           float *a, int *lda, float *beta, float *c, int *ldc);

int stbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           float *a, int *lda, float *x, int *incx);

int stbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           float *a, int *lda, float *x, int *incx);

int stpmv_(char *uplo, char *trans, char *diag, int *n, float *ap,
           float *x, int *incx);

int stpsv_(char *uplo, char *trans, char *diag, int *n, float *ap,
           float *x, int *incx);

int strmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, float *alpha, float *a, int *lda, float *b, 
           int *ldb);

int strmv_(char *uplo, char *trans, char *diag, int *n, float *a,
           int *lda, float *x, int *incx);

int strsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, float *alpha, float *a, int *lda, float *b, 
           int *ldb);

int strsv_(char *uplo, char *trans, char *diag, int *n, float *a,
           int *lda, float *x, int *incx);

int zaxpy_(int *n, dcomplex *ca, dcomplex *cx, int *incx, dcomplex *cy,
           int *incy);

int zcopy_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

int zdscal_(int *n, double *sa, dcomplex *cx, int *incx);

int zgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           dcomplex *alpha, dcomplex *a, int *lda, dcomplex *x, int *incx,
           dcomplex *beta, dcomplex *y, int *incy);

int zgemm_(char *transa, char *transb, int *m, int *n, int *k,
           dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, int *ldb,
           dcomplex *beta, dcomplex *c, int *ldc);

int zgemv_(char *trans, int *m, int *n, dcomplex *alpha, dcomplex *a,
           int *lda, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y,
           int *incy);

int zgerc_(int *m, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *a, int *lda);

int zgeru_(int *m, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *a, int *lda);

int zhbmv_(char *uplo, int *n, int *k, dcomplex *alpha, dcomplex *a,
           int *lda, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y,
           int *incy);

int zhemm_(char *side, char *uplo, int *m, int *n, dcomplex *alpha,
           dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *beta,
           dcomplex *c, int *ldc);

int zhemv_(char *uplo, int *n, dcomplex *alpha, dcomplex *a, int *lda,
           dcomplex *x, int *incx, dcomplex *beta, dcomplex *y, int *incy);

int zher_(char *uplo, int *n, double *alpha, dcomplex *x, int *incx,
          dcomplex *a, int *lda);

int zher2_(char *uplo, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *a, int *lda);

int zher2k_(char *uplo, char *trans, int *n, int *k, dcomplex *alpha,
            dcomplex *a, int *lda, dcomplex *b, int *ldb, double *beta,
            dcomplex *c, int *ldc);

int zherk_(char *uplo, char *trans, int *n, int *k, double *alpha,
           dcomplex *a, int *lda, double *beta, dcomplex *c, int *ldc);

int zhpmv_(char *uplo, int *n, dcomplex *alpha, dcomplex *ap, dcomplex *x,
           int *incx, dcomplex *beta, dcomplex *y, int *incy);

int zhpr_(char *uplo, int *n, double *alpha, dcomplex *x, int *incx,
          dcomplex *ap);

int zhpr2_(char *uplo, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *ap);

int zrotg_(dcomplex *ca, dcomplex *cb, double *c, dcomplex *s);

int zscal_(int *n, dcomplex *ca, dcomplex *cx, int *incx);

int zswap_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

int zsymm_(char *side, char *uplo, int *m, int *n, dcomplex *alpha,
           dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *beta,
           dcomplex *c, int *ldc);

int zsyr2k_(char *uplo, char *trans, int *n, int *k, dcomplex *alpha,
            dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *beta,
            dcomplex *c, int *ldc);

int zsyrk_(char *uplo, char *trans, int *n, int *k, dcomplex *alpha,
           dcomplex *a, int *lda, dcomplex *beta, dcomplex *c, int *ldc);

int ztbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           dcomplex *a, int *lda, dcomplex *x, int *incx);

int ztbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           dcomplex *a, int *lda, dcomplex *x, int *incx);

int ztpmv_(char *uplo, char *trans, char *diag, int *n, dcomplex *ap,
           dcomplex *x, int *incx);

int ztpsv_(char *uplo, char *trans, char *diag, int *n, dcomplex *ap,
           dcomplex *x, int *incx);

int ztrmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b,
           int *ldb);

int ztrmv_(char *uplo, char *trans, char *diag, int *n, dcomplex *a,
           int *lda, dcomplex *x, int *incx);

int ztrsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b,
           int *ldb);

int ztrsv_(char *uplo, char *trans, char *diag, int *n, dcomplex *a,
           int *lda, dcomplex *x, int *incx);
