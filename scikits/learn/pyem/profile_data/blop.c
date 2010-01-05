#include <math.h>
#include <stddef.h>

int compute(const double *in, size_t n, size_t d, const double* mu, double* out)
{
    size_t  i, j;
    double acc;

    for (i = 0; i < n; ++i) {
        acc = 0;
        for (j = 0; j < d; ++j) {
            acc += (in[i*d+j] - mu[j]) * (in[i*d+j] - mu[j]); 
        }
        out[i] = exp(acc);
    }

    return 0;
}

#if 0
int main(void) 
{
    const size_t n = 1e5;
    const size_t d = 30;
    size_t iter = 10, i;

    double  *in, *out;

    in = malloc(sizeof(*in) * n * d);
    out = malloc(sizeof(*out) * n * d);

    for (i = 0; i < iter; ++i) {
    }
    free(in);
    out(in);
}
#endif
