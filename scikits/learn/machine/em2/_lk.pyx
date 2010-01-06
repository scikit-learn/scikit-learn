import numpy as np
cimport numpy as c_np
cimport stdlib

cdef extern from "math.h":
    double exp(double)
    double log(double)

# XXX: find out why this one does not work for d > 1
#def quadform(c_np.ndarray x, c_np.ndarray mu, c_np.ndarray inva,
#             double fac, c_np.ndarray y):
#    """Compute (x-mu)'inva(x-mu) quadratic form for one set of K mean/inva
#    vectors of d dimension (inva is assumed to be diagonal, and inva contains
#    the diagonal items)."""
#    cdef double *raw_x, *raw_y, *raw_mu, *raw_inva
#    cdef int n, d, i
#    cdef double *cx
#
#    if not x.dtype == np.float64:
#        raise ValueError("Only float64 supported for now")
#
#    if not x.ndim == 2:
#        raise ValueError("Rank != 2 not supported yet")
#    n = x.shape[0]
#    d = x.shape[1]
#
#    if not inva.ndim == 1:
#        raise ValueError("inva has rank %d (expect 1)" % inva.ndim)
#    if not inva.shape[0] == d:
#        raise ValueError("inva has dim %d (expect %d)" % (inva.shape[0], d))
#    #if not fac.ndim == 1:
#    #    raise ValueError("fac has rank %d (expect 1)" % fac.ndim)
#    #if not fac.shape[0] == d:
#    #    raise ValueError("fac has dim %d (expect %d)" % (fac.shape[0], d))
#    if not mu.ndim == 1:
#        raise ValueError("mu has rank %d (expect 1)" % mu.ndim)
#    if not mu.shape[0] == d:
#        raise ValueError("mu has dim %d (expect %d)" % (mu.shape[0], d))
#
#    raw_x = <double*>x.data
#    raw_y = <double*>y.data
#    raw_mu = <double*>mu.data
#    raw_inva = <double*>inva.data
#
#    for i in range(n):
#        cx = raw_x + d * i
#        acc = 0
#        for j in range(d):
#            acc += raw_inva[j] * (cx[j] - raw_mu[j]) * (cx[j] - raw_mu[j])
#        y[i] = acc + fac
#
#    return y

def mquadform(c_np.ndarray x, c_np.ndarray mu, c_np.ndarray inva,
              c_np.ndarray fac, c_np.ndarray y):
    """Compute (x-mu)'inva(x-mu) quadratic form for a set of K mean/inva
    vectors of d dimension (inva is assumed to be diagonal, and inva contains
    the diagonal items)."""
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

def logsumexp(c_np.ndarray x, c_np.ndarray out):
    cdef double *raw_x, *raw_out
    cdef int n, d, no

    if not x.ndim == 2:
        raise ValueError("Rank != 2 not supported yet")

    if not x.dtype == np.float64:
        raise ValueError("Only float64 supported for now")

    n = x.shape[0]
    no = out.shape[0]
    d = x.shape[1]

    if not out.ndim == 1:
        raise ValueError("out should be a rank 1 array")
    if not no == n:
        raise ValueError(
            "out should have as many items as x's rows " \
            "(got %d, expected %d)" % (no, n))

    raw_x = <double*>x.data
    raw_out = <double*>out.data

    logsumexp_double(raw_x, n, d, raw_out)

cdef int logsumexp_double(double* x, int n, int d, double *out):
    cdef int i, j
    cdef double m, acc
    cdef double *cx, *cout

    for i in range(n):
        cx = x + i * d
        # For each row of x, compute m + log(e^{x[0]-m} + ... e^{x[d-1]-m})
        # where m is the maximum of x[0] ... x[d-1]
        m = cx[0]
        for j in range(1, d):
            if cx[j] > m:
                m = cx[j]
        acc = 0
        for j in range(d):
            acc += exp(cx[j] - m)
        out[i] = m + log(acc)

    return 0
