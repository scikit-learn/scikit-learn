#include <stddef.h>
#include <stdio.h>
#include <math.h>

/*
 * use va for caching, so its content is modified
 */
int gden_diag(const double *in, const size_t n, const size_t d, 
        const double* mu, double* va, double* out)
{
    size_t  j, i;
    double  fac = 1.0;
    double  acc, tmp;

    /*
     * Cache some precomputing
     */
    for(j = 0; j < d; ++j) {
        va[j]   = 0.5/va[j];
        fac     *= sqrt(va[j] / (M_PI) );
        va[j]   *= -1.0;
    }

    /*
     * actual computing
     */
    for(i = 0; i < n; ++i) {
        acc = 0;
        for(j = 0; j < d; ++j) {
            tmp = (in[i *d + j] - mu[j]);
            acc += tmp * tmp * va[j];
        }
        out[i]  = exp(acc) * fac;
    }
 
    return 0;
}
