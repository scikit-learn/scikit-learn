#include <math.h>
#include <cblas.h>

#ifdef _MSC_VER
# define inline __inline
#endif


/*
 * General Cholesky Delete.
 * Remove an element from the cholesky factorization
 * m = columns
 * n = rows
 *
 * TODO: put transpose as an option
 */
static inline void cholesky_delete_dbl(int m, int n, double *L, int go_out)
{
    double c, s;

    /* delete row go_out */
    double *L1 = L + (go_out * m);
    int i;
    for (i = go_out; i < n - 1; ++i) {
        cblas_dcopy (i + 2, L1 + m , 1, L1,  1);
        L1 += m;
    }

    L1 = L + (go_out * m);
    for (i=go_out; i < n - 1; ++i) {

        cblas_drotg(L1 + i, L1 + i + 1, &c, &s);
        if (L1[i] < 0) {
            /* Diagonals cannot be negative */
            L1[i] = fabs(L1[i]);
            c = -c;
            s = -s;
        }
        L1[i+1] = 0.; /* just for cleanup */
        L1 += m;

        cblas_drot(n - (i + 2), L1 + i, m, L1 + i + 1,
                   m, c, s);
    }
}


static inline void cholesky_delete_flt(int m, int n, float *L, int go_out)
{
    float c, s;

    /* delete row go_out */
    float *L1 = L + (go_out * m);
    int i;
    for (i = go_out; i < n - 1; ++i) {
        cblas_scopy (i + 2, L1 + m , 1, L1,  1);
        L1 += m;
    }

    L1 = L + (go_out * m);
    for (i=go_out; i < n - 1; ++i) {

        cblas_srotg(L1 + i, L1 + i + 1, &c, &s);
        if (L1[i] < 0) {
            /* Diagonals cannot be negative */
            L1[i] = fabsf(L1[i]);
            c = -c;
            s = -s;
        }
        L1[i+1] = 0.; /* just for cleanup */
        L1 += m;

        cblas_srot(n - (i + 2), L1 + i, m, L1 + i + 1,
                   m, c, s);
    }
}
