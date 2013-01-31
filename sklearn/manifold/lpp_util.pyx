'''
Created on 2012/11/30

@author: du
'''
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
def to_symmetrix_matrix(G, np.ndarray[int, ndim = 1] nzx, 
                   np.ndarray[int, ndim = 1] nzy):
    cdef int i, xx, yy, n_nz
    n_nz = len(nzx)
    for i in xrange(n_nz):
        xx = nzx[i]
        yy = nzy[i]
        if G[yy, xx] == 0:
            G[yy, xx] = G[xx, yy]
    return G