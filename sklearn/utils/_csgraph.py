"""Compressed Sparse graph algorithms"""
# Backported from scipy 0.9: scipy.sparse.csgraph

# Some compatibility fixes for scipy 0.6
# Fabian Pedregosa, October 2010

__docformat__ = "restructuredtext en"

__all__ = ['cs_graph_components']

import numpy as np

from .sparsetools import cs_graph_components as _cs_graph_components

from scipy.sparse import csr_matrix
from scipy.sparse import isspmatrix

_msg0 = 'x must be a symmetric square matrix!'
_msg1 = _msg0 + '(has shape %s)'


def cs_graph_components(x):
    """
    Determine connected compoments of a graph stored as a compressed sparse row
    or column matrix. For speed reasons, the symmetry of the matrix x is not
    checked.

    Parameters
    ----------
    x: ndarray-like, 2 dimensions, or sparse matrix
        The adjacency matrix of the graph. Only the upper triangular part
        is used.

    Returns
    -------
    n_components: int
        The number of connected components.
    label: ndarray (ints, 1 dimension):
        The label array of each connected component (-2 is used to
        indicate empty rows: 0 everywhere, including diagonal).

    Notes
    -----

    The matrix is assumed to be symmetric and the upper triangular part
    of the matrix is used. The matrix is converted to a CSR matrix unless
    it is already a CSR.

    Examples
    --------

    >>> from scipy.sparse import cs_graph_components
    >>> import numpy as np
    >>> D = np.eye(4)
    >>> D[0,1] = D[1,0] = 1
    >>> cs_graph_components(D)
    (3, array([0, 0, 1, 2]))
    >>> from scipy.sparse import dok_matrix
    >>> cs_graph_components(dok_matrix(D))
    (3, array([0, 0, 1, 2]))

    """
    try:
        shape = x.shape
    except AttributeError:
        raise ValueError(_msg0)

    if not ((len(x.shape) == 2) and (x.shape[0] == x.shape[1])):
        raise ValueError(_msg1 % x.shape)

    if isspmatrix(x):
        x = x.tocsr()
    else:
        x = csr_matrix(x)

    label = np.empty((shape[0],), dtype=x.indptr.dtype)

    n_components = _cs_graph_components(shape[0], x.indptr, x.indices, label)

    return n_components, label
