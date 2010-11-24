#include <math.h>
#include <cblas.h>

#ifdef _MSC_VER
    #define copysign _copysign
#endif

/* General Cholesky Delete.
 * Remove an element from the cholesky factorization
 * m = columns
 * n = rows
 *
 * TODO: put transpose as an option
 *
 */
int double_cholesky_delete (int m, int n, double *L, int go_out) {

    double c, s;

    /* delete row go_out */
    double * _L = L + (go_out * m);
    int i;
    for (i = go_out; i < n - 1; ++i) {
        cblas_dcopy (i + 2, _L + m , 1, _L,  1);
        _L += m;
    }

    _L = L + (go_out * m);
    for (i=go_out; i < n - 1; ++i) {

        cblas_drotg (_L + i, _L + i + 1, &c, &s);
        if (_L[i] < 0) {
            /* Diagonals cannot be negative */
            _L[i] = copysign(_L[i], 1.0);
            c = -c;
            s = -s;
        }
        _L[i+1] = 0.; /* just for cleanup */
        _L += m;

        cblas_drot (n - (i + 2), _L + i, m, _L + i + 1,
                    m, c, s);

    }

    return 0;
}


int float_cholesky_delete (int m, int n, float *L, int go_out) {

    float c, s;

    /* delete row go_out */
    float * _L = L + (go_out * m);
    int i;
    for (i = go_out; i < n - 1; ++i) {
        cblas_scopy (i + 2, _L + m , 1, _L,  1);
        _L += m;
    }

    _L = L + (go_out * m);
    for (i=go_out; i < n - 1; ++i) {

        cblas_srotg (_L + i, _L + i + 1, &c, &s);
        if (_L[i] < 0) {
            /* Diagonals cannot be negative */
            /* _L[i] = copysign(_L[i], 1.0); */
            c = -c;
            s = -s;
        }
        _L[i+1] = 0.; /* just for cleanup */
        _L += m;

        cblas_srot (n - (i + 2), _L + i, m, _L + i + 1,
                    m, c, s);

    }

    return 0;
}
