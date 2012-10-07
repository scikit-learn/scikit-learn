"""Algorithms for spectral clustering"""

# Author: Gael Varoquaux gael.varoquaux@normalesup.org
# License: BSD
import warnings

import numpy as np


from ..base import BaseEstimator, ClusterMixin
from ..utils import check_random_state
from ..utils.graph import graph_laplacian
from ..metrics.pairwise import rbf_kernel
from ..neighbors import kneighbors_graph
from .k_means_ import k_means

def _set_diag(laplacian, value):
    from scipy import sparse
    n_nodes = laplacian.shape[0]
    # We need to put the diagonal at zero
    if not sparse.isspmatrix(laplacian):
        laplacian.flat[::n_nodes + 1] = value
    else:
        laplacian = laplacian.tocoo()
        diag_idx = (laplacian.row == laplacian.col)
        laplacian.data[diag_idx] = value
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
    return laplacian


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

    mode: {None, 'arpack', 'lobpcg', or 'amg'}
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
    from scipy.sparse.linalg.eigen.lobpcg.lobpcg import symeig
    try:
        from pyamg import smoothed_aggregation_solver
    except ImportError:
        if mode == "amg":
            raise ValueError("The mode was set to 'amg', but pyamg is "
                             "not available.")

    random_state = check_random_state(random_state)

    n_nodes = adjacency.shape[0]
    # XXX: Should we check that the matrices given is symmetric
    if mode is None:
        mode = 'arpack'
    elif not mode in ('arpack', 'lobpcg', 'amg'):
        raise ValueError("Unknown value for mode: '%s'."
                         "Should be 'amg' or 'arpack'" % mode)
    laplacian, dd = graph_laplacian(adjacency,
                                    normed=True, return_diag=True)
    if (mode == 'arpack'
        or not sparse.isspmatrix(laplacian)
        or n_nodes < 5 * n_components):
        # lobpcg used with mode='amg' has bugs for low number of nodes
        laplacian = _set_diag(laplacian, 0)

        # Here we'll use shift-invert mode for fast eigenvalues
        # (see http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html
        #  for a short explanation of what this means)
        # Because the normalized Laplacian has eigenvalues between 0 and 2,
        # I - L has eigenvalues between -1 and 1.  ARPACK is most efficient
        # when finding eigenvalues of largest magnitude (keyword which='LM')
        # and when these eigenvalues are very large compared to the rest.
        # For very large, very sparse graphs, I - L can have many, many
        # eigenvalues very near 1.0.  This leads to slow convergence.  So
        # instead, we'll use ARPACK's shift-invert mode, asking for the
        # eigenvalues near 1.0.  This effectively spreads-out the spectrum
        # near 1.0 and leads to much faster convergence: potentially an
        # orders-of-magnitude speedup over simply using keyword which='LA'
        # in standard mode.
        try:
            lambdas, diffusion_map = eigsh(-laplacian, k=n_components,
                                        sigma=1.0, which='LM')
            embedding = diffusion_map.T[::-1] * dd
        except RuntimeError:
            # When submatrices are exactly singular, an LU decomposition
            # in arpack fails. We fallback to lobpcg
            mode = "lobpcg"

    if mode == 'amg':
        # Use AMG to get a preconditioner and speed up the eigenvalue
        # problem.
        laplacian = laplacian.astype(np.float)  # lobpcg needs native floats
        ml = smoothed_aggregation_solver(laplacian.tocsr())
        M = ml.aspreconditioner()
        X = random_state.rand(laplacian.shape[0], n_components)
        #X[:, 0] = 1. / dd.ravel()
        X[:, 0] = dd.ravel()
        lambdas, diffusion_map = lobpcg(laplacian, X, M=M, tol=1.e-12,
                                        largest=False)
        embedding = diffusion_map.T * dd
        if embedding.shape[0] == 1:
            raise ValueError
    elif mode == "lobpcg":
        laplacian = laplacian.astype(np.float)  # lobpcg needs native floats
        if n_nodes < 5 * n_components + 1:
            # lobpcg will fallback to symeig, so we short circuit it
            if sparse.isspmatrix(laplacian):
                laplacian = laplacian.todense()
            lambdas, diffusion_map = symeig(laplacian)
            embedding = diffusion_map.T[:n_components] * dd
        else:
            laplacian = laplacian.astype(np.float)  # lobpcg needs native floats
            laplacian = _set_diag(laplacian, 1)
            # We increase the number of eigenvectors requested, as lobpcg
            # doesn't behave well in low dimension
            X = random_state.rand(laplacian.shape[0], n_components + 1)
            X[:, 0] = dd.ravel()
            lambdas, diffusion_map = lobpcg(laplacian, X, tol=1e-15,
                                            largest=False, maxiter=2000)
            embedding = diffusion_map.T[:n_components] * dd
            if embedding.shape[0] == 1:
                raise ValueError
    return embedding


def spectral_clustering(affinity, n_clusters=8, n_components=None, mode=None,
                        random_state=None, n_init=10, k=None):
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

    n_clusters: integer, optional
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

    References
    ----------

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
    if not k is None:
        warnings.warn("'k' was renamed to n_clusters", DeprecationWarning)
        n_clusters = k
    random_state = check_random_state(random_state)
    n_components = n_clusters if n_components is None else n_components
    maps = spectral_embedding(affinity, n_components=n_components,
                              mode=mode, random_state=random_state)
    maps = maps[1:]
    _, labels, _ = k_means(maps.T, n_clusters, random_state=random_state,
                    n_init=n_init)
    return labels


class SpectralClustering(BaseEstimator, ClusterMixin):
    """Apply k-means to a projection to the normalized laplacian

    In practice Spectral Clustering is very useful when the structure of
    the individual clusters is highly non-convex or more generally when
    a measure of the center and spread of the cluster is not a suitable
    description of the complete cluster. For instance when clusters are
    nested circles on the 2D plan.

    If affinity is the adjacency matrix of a graph, this method can be
    used to find normalized graph cuts.

    When calling ``fit``, an affinity matrix is constructed using either the
    Gaussian (aka RBF) kernel of the euclidean distanced ``d(X, X)``::

            np.exp(-gamma * d(X,X) ** 2)

    or a k-nearest neighbors connectivity matrix.

    Alternatively, using ``precomputed``, a user-provided affinity
    matrix can be used.

    Parameters
    -----------
    n_clusters : integer, optional
        The dimension of the projection subspace.

    affinity: string, 'nearest_neighbors', 'rbf' or 'precomputed'

    gamma: float
        Scaling factor of Gaussian (rbf) affinity kernel. Ignored for
        ``affinity='nearest_neighbors'``.

    n_neighbors: integer
        Number of neighbors to use when constructing the affinity matrix using
        the nearest neighbors method. Ignored for ``affinity='rbf'``.

    mode: {None, 'arpack' or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities

    random_state : int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization
        of the lobpcg eigen vectors decomposition when mode == 'amg'
        and by the K-Means initialization.

    n_init : int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.


    Attributes
    ----------
    `affinity_matrix_` : array-like, shape (n_samples, n_samples)
        Affinity matrix used for clustering. Available only if after calling
        ``fit``.

    `labels_` :
        Labels of each point

    Notes
    -----
    If you have an affinity matrix, such as a distance matrix,
    for which 0 means identical elements, and high values means
    very dissimilar elements, it can be transformed in a
    similarity matrix that is well suited for the algorithm by
    applying the Gaussian (RBF, heat) kernel::

        np.exp(- X ** 2 / (2. * delta ** 2))

    Another alternative is to take a symmetric version of the k
    nearest neighbors connectivity matrix of the points.

    If the pyamg package is installed, it is used: this greatly
    speeds up computation.

    References
    ----------

    - Normalized cuts and image segmentation, 2000
      Jianbo Shi, Jitendra Malik
      http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.160.2324

    - A Tutorial on Spectral Clustering, 2007
      Ulrike von Luxburg
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.165.9323
    """

    def __init__(self, n_clusters=8, mode=None, random_state=None, n_init=10,
            gamma=1., affinity='rbf', n_neighbors=10, k=None,
            precomputed=False):
        if not k is None:
            warnings.warn("'k' was renamed to n_clusters", DeprecationWarning)
            n_clusters = k
        self.n_clusters = n_clusters
        self.mode = mode
        self.random_state = random_state
        self.n_init = n_init
        self.gamma = gamma
        self.affinity = affinity
        self.n_neighbors = n_neighbors

    def fit(self, X):
        """Creates an affinity matrix for X using the selected affinity,
        then applies spectral clustering to this affinity matrix.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            OR, if affinity==`precomputed`, a precomputed affinity
            matrix of shape (n_samples, n_samples)
        """
        if X.shape[0] == X.shape[1] and self.affinity != "precomputed":
            warnings.warn("The spectral clustering API has changed. ``fit``"
                    "now constructs an affinity matrix from data. To use "
                    "a custom affinity matrix, set ``affinity=precomputed``.")

        if self.affinity == 'rbf':
            self.affinity_matrix_ = rbf_kernel(X, gamma=self.gamma)

        elif self.affinity == 'nearest_neighbors':
            connectivity = kneighbors_graph(X, n_neighbors=self.n_neighbors)
            self.affinity_matrix_ = 0.5 * (connectivity + connectivity.T)
        elif self.affinity == 'precomputed':
            self.affinity_matrix_ = X
        else:
            raise ValueError("Invalid 'affinity'. Expected 'rbf', "
                "'nearest_neighbors' or 'precomputed', got '%s'."
                % self.affinity_matrix)

        self.random_state = check_random_state(self.random_state)
        self.labels_ = spectral_clustering(self.affinity_matrix_,
                n_clusters=self.n_clusters, mode=self.mode,
                random_state=self.random_state, n_init=self.n_init)
        return self

    @property
    def _pairwise(self):
        return self.affinity == "precomputed"
