"""Power Iteration Clustering

Scalable alternative to Spectral Clustering for small number of centers.
"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD

import numpy as np

from ..utils.extmath import safe_sparse_dot


def power_iteration_clustering(affinity, k=8, n_vectors=1, tol=1e-5,
                               normalize=True):
    """Power Iteration Clustering: simple variant of spectral clustering

    One or more random vectors are multiplied several times to the
    row normalized affinity matrix so as to reach a local convergence
    (early stopping before reaching the convergence to the first
    eigen-vector).

    This process imprints the features of the major eigen-vectors into
    the vectors to make them suitable as clustering features for a couple
    of fast K-Means.

    Parameters
    ----------

      TODO

    Returns
    --------
    labels: array of integer, shape: (n_samples, k)
        The cluster label assignement for each sample.


    Reference
    ---------

    W. Cohen, F. Lin, Power Iteration Clustering, ICML 2010
    http://www.cs.cmu.edu/~wcohen/postscript/icml2010-pic-final.pdf

    Complexity
    ----------

    TODO

    """

    if normalize:
        # row normalize the affinity matrix
        sums = affinity.sum(axis=1)
        nnzeros = np.where(scales > 0)
        scales[nnzeros] = 1 / sums[nnzeros]

        if hasattr(affinity, 'tocsr'):
            # inplace row normalization for sparse matrices

            # late import of scipy.sparse for performance
            from scipy.sparse.sparsetools import csr_scale_rows

            # ensure the sparse matrix is in Compressed Sparse Rows format
            affinity = affinity.tocsr()

            # convert matrix to array
            scales = np.array(scales)[0]

            # inplace rescaling of the CSR matrix
            csr_scale_rows(affinity.shape[0], affinity.shape[1],
                           affinity.indptr, affinity.indices,
                           affinity.data, scales)
        else:
            # inplace row normalization for ndarray
            affinity = np.asanyarray(affinity)
            affinity /= scales.reshape((s.shape[0], -1))

    #TODO




