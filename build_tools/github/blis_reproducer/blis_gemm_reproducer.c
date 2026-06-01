#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum CBLAS_ORDER { CblasRowMajor = 101, CblasColMajor = 102 };
enum CBLAS_TRANSPOSE { CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113 };

extern void cblas_dgemm(
    enum CBLAS_ORDER order,
    enum CBLAS_TRANSPOSE transa,
    enum CBLAS_TRANSPOSE transb,
    int m,
    int n,
    int k,
    double alpha,
    const double *a,
    int lda,
    const double *b,
    int ldb,
    double beta,
    double *c,
    int ldc
);

static double deterministic_value(int i, int j) {
    return sin((double)(i + 1) * 0.173) + cos((double)(j + 1) * 0.113);
}

int main(void) {
    const int m = 192;
    const int n = 160;
    const int k = 128;
    const double alpha = 1.0;
    const double beta = 0.0;
    const double atol = 1e-11;
    const double rtol = 1e-11;

    double *a = calloc((size_t)m * (size_t)k, sizeof(double));
    double *b = calloc((size_t)k * (size_t)n, sizeof(double));
    double *c = calloc((size_t)m * (size_t)n, sizeof(double));
    long double *c_ref = calloc((size_t)m * (size_t)n, sizeof(long double));

    if (!a || !b || !c || !c_ref) {
        fprintf(stderr, "allocation failure\n");
        free(a);
        free(b);
        free(c);
        free(c_ref);
        return 2;
    }

    for (int i = 0; i < m; i++) {
        for (int p = 0; p < k; p++) {
            a[(size_t)i * (size_t)k + (size_t)p] = deterministic_value(i, p);
        }
    }
    for (int p = 0; p < k; p++) {
        for (int j = 0; j < n; j++) {
            b[(size_t)p * (size_t)n + (size_t)j] = deterministic_value(p + 7, j + 11);
        }
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, k, b, n, beta, c, n);

    long double max_abs_err = 0.0L;
    long double max_rel_err = 0.0L;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            long double acc = 0.0L;
            for (int p = 0; p < k; p++) {
                acc += (long double)a[(size_t)i * (size_t)k + (size_t)p] *
                       (long double)b[(size_t)p * (size_t)n + (size_t)j];
            }
            c_ref[(size_t)i * (size_t)n + (size_t)j] = acc;
            long double diff = fabsl((long double)c[(size_t)i * (size_t)n + (size_t)j] - acc);
            if (diff > max_abs_err) {
                max_abs_err = diff;
            }
            long double denom = fabsl(acc);
            if (denom < 1e-15L) {
                denom = 1e-15L;
            }
            long double rel = diff / denom;
            if (rel > max_rel_err) {
                max_rel_err = rel;
            }
        }
    }

    printf("max_abs_err=%.6Le\n", max_abs_err);
    printf("max_rel_err=%.6Le\n", max_rel_err);

    free(a);
    free(b);
    free(c);
    free(c_ref);

    if (max_abs_err > atol && max_rel_err > rtol) {
        fprintf(stderr, "FAIL: BLIS GEMM mismatch above tolerance\n");
        return 1;
    }

    printf("PASS\n");
    return 0;
}
