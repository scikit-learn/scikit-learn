"""Algorithms for spectral clustering"""

# Author: Gael Varoquaux gael.varoquaux@normalesup.org
# License: BSD
import warnings

import numpy as np


from ..base import BaseEstimator
from ..utils import check_random_state
from ..utils.graph import graph_laplacian
from .k_means_ import k_means


def spectral_embedding(adjacency, n_components=8, mode=None,
                       random_state=None):
    """Project the sample on the first eigen vectors of the graph Laplacian

    The adjacency matrix is used to compute a normalized graph Laplacian
    whose spectrum (especially the eigen vectors associated to the
    smallest eigen values) has an interpretation in terms of minimal
    number of cuts necessary to split the graph into comparably sized
    components.

    This embedding can also 'work' even if the ``adjacency`` variable is
    not strictly the adjacency matrix of a graph but more generally
    an affinity or similarity matrix between samples (for instance the
    heat kernel of a euclidean distance matrix or a k-NN matrix).

    However care must taken to always make the affinity matrix symmetric
    so that the eigen vector decomposition works as expected.

    Parameters
    -----------
    adjacency: array-like or sparse matrix, shape: (n_samples, n_samples)
        The adjacency matrix of the graph to embed.

    n_components: integer, optional
        The dimension of the projection subspace.

    mode: {None, 'arpack' or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities

    random_state: int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization of the
        lobpcg eigen vectors decomposition when mode == 'amg'. By default
        arpack is used.

    Returns
    --------
    embedding: array, shape: (n_samples, n_components)
        The reduced samples

    Notes
    ------
    The graph should contain only one connected component, elsewhere the
    results make little sense.
    """

    from scipy import sparse
    from ..utils.arpack import eigsh
    from scipy.sparse.linalg import lobpcg
    try:
        from pyamg import smoothed_aggregation_solver
        amg_loaded = True
    except ImportError:
        amg_loaded = False

    random_state = check_random_state(random_state)

    n_nodes = adjacency.shape[0]
    # XXX: Should we check that the matrices given is symmetric
    if not amg_loaded:
        warnings.warn('pyamg not available, using scipy.sparse')
    if mode is None:
        mode = 'arpack'
    laplacian, dd = graph_laplacian(adjacency,
                                    normed=True, return_diag=True)
    if (mode == 'arpack'
        or not sparse.isspmatrix(laplacian)
        or n_nodes < 5 * n_components):
        # lobpcg used with mode='amg' has bugs for low number of nodes

        # We need to put the diagonal at zero
        if not sparse.isspmatrix(laplacian):
            laplacian[::n_nodes + 1] = 0
        else:
            laplacian = laplacian.tocoo()
            diag_idx = (laplacian.row == laplacian.col)
            laplacian.data[diag_idx] = 0
            # If the matrix has a small number of diagonals (as in the
            # case of structured matrices comming from images), the
            # dia format might be best suited for matvec products:
            n_diags = np.unique(laplacian.row - laplacian.col).size
            if n_diags <= 7:
                # 3 or less outer diagonals on each side
                laplacian = laplacian.todia()
            else:
                # csr has the fastest matvec and is thus best suited to
                # arpack
                laplacian = laplacian.tocsr()
        lambdas, diffusion_map = eigsh(-laplacian, k=n_components,
                                        which='LA')
        embedding = diffusion_map.T[::-1] * dd
    elif mode == 'amg':
        # Use AMG to get a preconditioner and speed up the eigenvalue
        # problem.
        laplacian = laplacian.astype(np.float)  # lobpcg needs native floats
        ml = smoothed_aggregation_solver(laplacian.tocsr())
        X = random_state.rand(laplacian.shape[0], n_components)
        X[:, 0] = 1. / dd.ravel()
        M = ml.aspreconditioner()
        lambdas, diffusion_map = lobpcg(laplacian, X, M=M, tol=1.e-12,
                                        largest=False)
        embedding = diffusion_map.T * dd
        if embedding.shape[0] == 1:
            raise ValueError
    else:
        raise ValueError("Unknown value for mode: '%s'."
                         "Should be 'amg' or 'arpack'" % mode)
    return embedding


def spectral_clustering(affinity, k=8, n_components=None, mode=None,
                        random_state=None, n_init=10):
    """Apply k-means to a projection to the normalized laplacian

    In practice Spectral Clustering is very useful when the structure of
    the individual clusters is highly non-convex or more generally when
    a measure of the center and spread of the cluster is not a suitable
    description of the complete cluster. For instance when clusters are
    nested circles on the 2D plan.

    If affinity is the adjacency matrix of a graph, this method can be
    used to find normalized graph cuts.

    Parameters
    -----------
    affinity: array-like or sparse matrix, shape: (n_samples, n_samples)
        The affinity matrix describing the relationship of the samples to
        embed. **Must be symetric**.

        Possible examples:
          - adjacency matrix of a graph,
          - heat kernel of the pairwise distance matrix of the samples,
          - symmetic k-nearest neighbours connectivity matrix of the samples.

    k: integer, optional
        Number of clusters to extract.

    n_components: integer, optional, default is k
        Number of eigen vectors to use for the spectral embedding

    mode: {None, 'arpack' or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities

    random_state: int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization
        of the lobpcg eigen vectors decomposition when mode == 'amg'
        and by the K-Means initialization.

    n_init: int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    Returns
    -------
    labels: array of integers, shape: n_samples
        The labels of the clusters.

    centers: array of integers, shape: k
        The indices of the cluster centers

    Notes
    -----
    **References**:

    - Normalized cuts and image segmentation, 2000
      Jianbo Shi, Jitendra Malik
      http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.160.2324

    - A Tutorial on Spectral Clustering, 2007
      Ulrike von Luxburg
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.165.9323

    Notes
    ------
    The graph should contain only one connect component, elsewhere
    the results make little sense.

    This algorithm solves the normalized cut for k=2: it is a
    normalized spectral clustering.
    """
    random_state = check_random_state(random_state)
    n_components = k if n_components is None else n_components
    maps = spectral_embedding(affinity, n_components=n_components,
                              mode=mode, random_state=random_state)
    maps = maps[1:]
    _, labels, _ = k_means(maps.T, k, random_state=random_state,
                    n_init=n_init)
    return labels


class SpectralClustering(BaseEstimator):
    """Apply k-means to a projection to the normalized laplacian

    In practice Spectral Clustering is very useful when the structure of
    the individual clusters is highly non-convex or more generally when
    a measure of the center and spread of the cluster is not a suitable
    description of the complete cluster. For instance when clusters are
    nested circles on the 2D plan.

    If affinity is the adjacency matrix of a graph, this method can be
    used to find normalized graph cuts.

    Parameters
    -----------
    k: integer, optional
        The dimension of the projection subspace.

    mode: {None, 'arpack' or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities

    random_state: int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization
        of the lobpcg eigen vectors decomposition when mode == 'amg'
        and by the K-Means initialization.

    n_init: int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    Attributes
    ----------

    `labels_` :
        Labels of each point

    Notes
    -----
    **References**:

    - Normalized cuts and image segmentation, 2000
      Jianbo Shi, Jitendra Malik
      http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.160.2324

    - A Tutorial on Spectral Clustering, 2007
      Ulrike von Luxburg
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.165.9323
    """

    def __init__(self, k=8, mode=None, random_state=None, n_init=10):
        self.k = k
        self.mode = mode
        self.random_state = random_state
        self.n_init = n_init

    def fit(self, X):
        """Compute the spectral clustering from the affinity matrix

        Parameters
        -----------
        X: array-like or sparse matrix, shape: (n_samples, n_samples)
            An affinity matrix describing the pairwise similarity of the
            data. If can also be an adjacency matrix of the graph to embed.
            X must be symmetric and its entries must be positive or
            zero. Zero means that elements have nothing in common,
            whereas high values mean that elements are strongly similar.

        Notes
        ------
        If you have an affinity matrix, such as a distance matrix,
        for which 0 means identical elements, and high values means
        very dissimilar elements, it can be transformed in a
        similarity matrix that is well suited for the algorithm by
        applying the gaussian (heat) kernel::

            np.exp(- X ** 2 / (2. * delta ** 2))

        Another alternative is to take a symmetric version of the k
        nearest neighbors connectivity matrix of the points.

        If the pyamg package is installed, it is used: this greatly
        speeds up computation.
        """
        self.random_state = check_random_state(self.random_state)
        self.labels_ = spectral_clustering(X, k=self.k, mode=self.mode,
                                           random_state=self.random_state,
                                           n_init=self.n_init)
        return self
