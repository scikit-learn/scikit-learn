import numpy as np
cimport numpy as c_np
cimport stdlib

def quadform(c_np.ndarray x, c_np.ndarray mu, c_np.ndarray inva, 
             c_np.ndarray fac, c_np.ndarray y):
    cdef double *raw_x, *raw_y, *raw_mu, *raw_inva, *raw_fac
    cdef int n, k, d

    if not x.dtype == np.float64:
        raise ValueError("Only float64 supported for now")

    if not x.ndim == 2:
        raise ValueError("Rank != 2 not supported yet")
    n = x.shape[0]
    d = x.shape[1]
    k = mu.shape[0]

    raw_x = <double*>x.data
    raw_y = <double*>y.data
    raw_mu = <double*>mu.data
    raw_inva = <double*>inva.data
    raw_fac = <double*>fac.data

    quadform_double(raw_x, n, d, raw_mu, raw_inva, raw_fac, k, raw_y)
    return y

# y values are assumed to be set to 0
cdef int quadform_double(double* x, int n, int d, double *mu, double* inva, 
                         double *fac, int k, double*y):
    cdef int i, j, c
    cdef double acc, *cx, *cmu, *cy, *cinva

    for i in range(n):
        cx = x + d * i
        cy = y + k * i
        for c in range(k):
            cmu = mu + d * c
            cinva = inva + d * c
            acc = 0
            for j in range(d):
                acc += cinva[j] * (cx[j] - cmu[j]) * (cx[j] - cmu[j])
            cy[c] = acc + fac[c]

    return 0
