# Last Change: Thu Jul 13 02:00 PM 2006 J

cimport c_numpy
cimport c_python
import numpy

c_numpy.import_array()

# pyrex version of _vq. Much faster in high dimension/high number of K 
# (ie K > 3-4)
def _vq(data, init):
    (n, d)  = data.shape
    label   = numpy.zeros(n, int)
    _imp_vq(data, init, label)

    return label

def _imp_vq(c_numpy.ndarray data, c_numpy.ndarray init, c_numpy.ndarray label):
    cdef int n
    cdef int d
    cdef int nc
    cdef int i
    cdef int j
    cdef int k
    cdef int *labeld
    cdef double *da, *code
    cdef double dist
    cdef double acc

    n   = data.dimensions[0]
    d   = data.dimensions[1]
    nc  = init.dimensions[0]

    if not data.dtype == numpy.dtype(numpy.float64):
        print '_vq not (yet) implemented for dtype %s'%dtype.name
        return
    da  = (<double*>data.data)

    if not init.dtype == numpy.dtype(numpy.float64):
        print '_vq not (yet) implemented for dtype %s'%dtype.name
        return
    code    = (<double*>init.data)

    if not label.dtype == numpy.dtype(numpy.int32):
        print '_vq not (yet) implemented for dtype %s'%dtype.name
        return
    labeld  = (<int*>label.data)

    for i from 0<=i<n:
        acc = 0
        lab = 0
        for j from 0<=j<d:
            acc = acc + (da[i * d + j] - code[j]) * \
                (da[i * d + j] - code[j])
        dist    = acc
        for k from 1<=k<nc:
            acc     = 0
            for j from 0<=j<d:
                acc = acc + (da[i * d + j] - code[k * d + j]) * \
                    (da[i * d + j] - code[k * d + j])
            if acc < dist:
                dist    = acc
                lab     = k
        labeld[i]   = lab

    return lab

# # Scalar gaussian a posteriori pdf
# def _scalar_gauss_den(c_numpy.ndarray x, double mu, 
#         double va, int log):
#     """ This function is the actual implementation
#     of gaussian pdf in scalar case. It assumes all args
#     are conformant, so it should not be used directly
#     
#     ** Expect centered data (ie with mean removed) **
# 
#     Call gauss_den instead"""
#     cdef int d, n
#     cdef double inva, fac
#     cdef double* data
# 
#     if not x.dtype == numpy.dtype(numpy.float64):
#         print '_scalar_gauss_den not (yet) implemented for dtype %s'%dtype.name
#         return
#     data  = (<double*>x.data)
# 
#     d       = x.dimensions[1]
#     n       = x.dimensions[0]
#     inva    = 1/va
#     fac     = (2*N.pi) ** (-d/2.0) * N.sqrt(inva)
#     y       = (x ** 2) * -0.5 * inva
#     if not log:
#         y   = fac * N.exp(y)
#     else:
#         y   = y + log(fac)
# 
#     return y
#     
