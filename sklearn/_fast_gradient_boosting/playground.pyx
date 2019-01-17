import numpy as np
from cython.parallel import prange


def wrapper():
    print('in')
    a = np.random.uniform(0, 100, size=(100, 100)).astype(np.int32)
    g(a)

cdef int f(int [:] a) nogil:
    return 3

cdef int g(int [:, :] a) nogil:

    cdef:
        int i

    for i in range(a.shape[0]):
        f(a[i])