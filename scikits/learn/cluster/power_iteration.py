"""Power Iteration Clustering

Scalable alternative to Spectral Clustering for small number of centers.
"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD

import numpy as np

from .k_means_ import k_means
from ..utils.extmath import safe_sparse_dot


def power_iteration_clustering(affinity, k=8, n_vectors=1, tol=1e-5,
                               rng=0, max_iter=1000, verbose=False):
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
    if rng is None:
        rng = np.random.RandomState()
    elif isinstance(rng, int):
        rng = np.random.RandomState(rng)

    if not hasattr(affinity, 'todense'):
        # this is not a sparse matrix: check that this is an array like
        affinity = np.asanyarray(affinity)

    # row normalize the affinity matrix
    sums = affinity.sum(axis=1)
    volume = sums.sum()

    scales = sums.copy()
    nnzeros = np.where(scales > 0)
    scales[nnzeros] = 1 / scales[nnzeros]

    if hasattr(affinity, 'tocsr'):
        # inplace row normalization for sparse matrices

        # late import of scipy.sparse for performance
        from scipy.sparse.sparsetools import csr_scale_rows

        # ensure the sparse matrix is in Compressed Sparse Rows format
        normalized = affinity.tocsr()
        if normalized is affinity:
            normalized = normalized.copy()

        # convert matrices to arrays
        scales = np.array(scales)[0]
        sums = np.array(sums)[0]

        # inplace rescaling of the CSR matrix
        csr_scale_rows(normalized.shape[0], normalized.shape[1],
                       normalized.indptr, normalized.indices,
                       normalized.data, scales)
    else:
        # inplace row normalization for ndarray
        normalized = affinity / scales[:, np.newaxis]

    n_samples = affinity.shape[0]

    if n_vectors == 1:
        # initialize a single vector deterministically
        vectors = (sums / volume).reshape((n_vectors, n_samples))
    else:
        # random init
        vectors = rng.normal(size=(n_vectors, n_samples))

    previous_vectors = vectors.copy()
    delta = 500

    for i in range(max_iter):

        previous_vectors[:] = vectors
        previous_delta = delta

        vectors[:] = safe_sparse_dot(normalized, vectors.T).T
        vectors /= np.abs(vectors).sum(axis=1)[:, np.newaxis]

        delta = np.abs(previous_vectors - vectors).mean()

        if verbose and i % 10 == 0:
            print "Power Iteration %04d/%04d: delta=%f" % (
                i + 1, max_iter, delta)

        if np.abs(previous_delta - delta) < tol:
            break
    if verbose:
        print "Converged at iteration: %04d/%04d with delta=%f" % (
            i + 1, max_iter, delta)

    # TODO: pass the rng correctly
    _, labels, inertia = k_means(vectors.T, k, verbose=verbose)
    return labels, inertia, vectors

