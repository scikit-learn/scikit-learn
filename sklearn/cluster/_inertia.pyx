import numpy as np
cimport numpy as np
cimport cython
ctypedef np.float64_t DOUBLE
ctypedef np.int_t INT


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)

def compute_ward_dist(np.ndarray[DOUBLE, ndim=1] m_1,\
                    np.ndarray[DOUBLE, ndim=2] m_2,\
                    np.ndarray[INT, ndim=1] coord_row,
                    np.ndarray[INT, ndim=1] coord_col,\
                    np.ndarray[DOUBLE, ndim=1] res):
    """ Compute distance of two cluster according to Ward's criterion.
    
    According to Ward's criterion, at each step during agglomerative clustering,
    the pair of clusters with minimum cluster distance are merged, where the
    cluster distance is the increase of within-cluster variance by merging
    two clusters.
    
    This increase in within-cluster variance when merging clusters A and B
    can be computed efficiently in constant time by using the formula 
    d(A,B) = d(c_a,c_b)**2/(1/n_a + 1/n_b)
    where c_a denotes the centroid of cluster A and n_a the number of elements
    in cluster A.
    
    Reference: https://secure.wikimedia.org/wikipedia/en/wiki/Nearest-neighbor_chain_algorithm
    """
    cdef int size_max = coord_row.shape[0]
    cdef int n_features = m_2.shape[1]
    cdef int i, j, row, col
    cdef DOUBLE pa, n
    for i in range(size_max):
        row = coord_row[i]
        col = coord_col[i]
        n = (m_1[row] * m_1[col]) / (m_1[row] + m_1[col])
        pa = 0.
        for j in range(n_features):
            pa += (m_2[row, j] / m_1[row] - m_2[col, j] / m_1[col])**2
        res[i] = pa * n
    return res
