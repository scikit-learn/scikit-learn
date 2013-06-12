import numpy as np
from collections import namedtuple

SparseBicluster = namedtuple('SparseBicluster', ('rows', 'columns'))

def to_sparse(rows, cols):
    """Transforms dense representation to sparse representation."""
    if rows.shape[0] != cols.shape[0]:
        raise Exception('number of biclusters is inconsistent')
    n_clusters = rows.shape[0]

    result = []
    for c in range(n_clusters):
        result.append(SparseBicluster(np.flatnonzero(rows[c]), np.flatnonzero(cols[c])))
    return result
