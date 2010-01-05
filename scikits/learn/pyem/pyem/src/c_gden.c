#include <stddef.h>
#include <stdio.h>
#include <math.h>

int gden_diag(const double *in, const size_t n, const size_t d, 
        const double* mu, const double* inva, double* out)
{
    size_t  j, i;
    double  fac = 1.0;
    double  acc, tmp;

    /*
     * Cache some precomputing
     */
    for(j = 0; j < d; ++j) {
        fac     *= sqrt(inva[j] / (2 * M_PI) );
        //printf("inva[%d] is %f, fac %f\n", j, inva[j], fac);
    }

    /*
     * actual computing
     */
    for(i = 0; i < n; ++i) {
        acc = 0;
        for(j = 0; j < d; ++j) {
            tmp = (in[d * i + j] - mu[j]);
            acc += tmp * tmp * -0.5 * inva[j];
        }
        out[i]  = exp(acc) * fac;
    }

    return 0;
}
