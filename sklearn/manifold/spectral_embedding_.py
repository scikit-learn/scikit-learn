"""Spectral Embedding"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Wei LI <kuantkid@gmail.com>
# License: BSD 3 clause

import warnings
import numpy as np

from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import lobpcg

from ..base import BaseEstimator
from ..externals import six
from ..utils import check_random_state
from ..utils.validation import check_array
from ..utils.graph import graph_laplacian
from ..utils.sparsetools import connected_components
from ..utils.arpack import eigsh
from ..metrics.pairwise import rbf_kernel
from ..neighbors import kneighbors_graph


def _graph_connected_component(graph, node_id):
    """Find the largest graph connected components the contains one
    given node

    Parameters
    ----------
    graph : array-like, shape: (n_samples, n_samples)
        adjacency matrix of the graph, non-zero weight means an edge
        between the nodes

    node_id : int
        The index of the query node of the graph

    Returns
    -------
    connected_components_matrix : array-like, shape: (n_samples,)
        An array of bool value indicates the indexes of the nodes
        belong to the largest connected components of the given query
        node
    """
    connected_components_matrix = np.zeros(shape=(graph.shape[0]), dtype=np.bool)
    connected_components_matrix[node_id] = True
    n_node = graph.shape[0]
    for i in range(n_node):
        last_num_component = connected_components_matrix.sum()
        _, node_to_add = np.where(graph[connected_components_matrix] != 0)
        connected_components_matrix[node_to_add] = True
        if last_num_component >= connected_components_matrix.sum():
            break
    return connected_components_matrix


def _graph_is_connected(graph):
    """ Return whether the graph is connected (True) or Not (False)

    Parameters
    ----------
    graph : array-like or sparse matrix, shape: (n_samples, n_samples)
        adjacency matrix of the graph, non-zero weight means an edge
        between the nodes

    Returns
    -------
    is_connected : bool
        True means the graph is fully connected and False means not
    """
    if sparse.isspmatrix(graph):
        # sparse graph, find all the connected components
        n_connected_components, _ = connected_components(graph)
        return n_connected_components == 1
    else:
        # dense graph, find all connected components start from node 0
        return _graph_connected_component(graph, 0).sum() == graph.shape[0]


def _set_diag(laplacian, value):
    """Set the diagonal of the laplacian matrix and convert it to a
    sparse format well suited for eigenvalue decomposition

    Parameters
    ----------
    laplacian : array or sparse matrix
        The graph laplacian
    value : float
        The value of the diagonal

    Returns
    -------
    laplacian : array or sparse matrix
        An array of matrix in a form that is well suited to fast
        eigenvalue decomposition, depending on the band width of the
        matrix.
    """
    n_nodes = laplacian.shape[0]
    # We need all entries in the diagonal to values
    if not sparse.isspmatrix(laplacian):
        laplacian.flat[::n_nodes + 1] = value
    else:
        laplacian = laplacian.tocoo()
        diag_idx = (laplacian.row == laplacian.col)
        laplacian.data[diag_idx] = value
        # If the matrix has a small number of diagonals (as in the
        # case of structured matrices coming from images), the
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


def _solve_eigenvalue_problem(adjacency, n_components=1, eigen_solver=None,
                              random_state=None, eigen_tol=0.0,
                              norm_laplacian=True, mode=None):
    """Project the sample on the first eigen vectors of the graph Laplacian.

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
    ----------
    adjacency : array-like or sparse matrix, shape: (n_samples, n_samples)
        The adjacency matrix of the graph to embed.

    n_components : integer, optional
        The dimension of the projection subspace.

    eigen_solver : {None, 'arpack', 'lobpcg', or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities.

    random_state : int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization of the
        lobpcg eigen vectors decomposition when eigen_solver == 'amg'.
        By default, arpack is used.

    eigen_tol : float, optional, default=0.0
        Stopping criterion for eigendecomposition of the Laplacian matrix
        when using arpack eigen_solver.

    drop_first : bool, optional, default=True
        Whether to drop the first eigenvector. For spectral embedding, this
        should be True as the first eigenvector should be constant vector for
        connected graph, but for spectral clustering, this should be kept as
        False to retain the first eigenvector.

    Returns
    -------
    lambdas : array, shape=(n_components,)
        The eigenvalues of the Laplacian
    vectors : array, shape=(n_components, n_samples)
        The eigenvectors of the Laplacian
    degrees : array, shape=(n_samples,)
        The degrees of the graph

    Notes
    -----
    Spectral embedding is most useful when the graph has one connected
    component. If there graph has many components, the first few eigenvectors
    will simply uncover the connected components of the graph.

    References
    ----------
    * http://en.wikipedia.org/wiki/LOBPCG

    * Toward the Optimal Preconditioned Eigensolver: Locally Optimal
      Block Preconditioned Conjugate Gradient Method
      Andrew V. Knyazev
      http://dx.doi.org/10.1137%2FS1064827500366124
    """

    try:
        from pyamg import smoothed_aggregation_solver
    except ImportError:
        if eigen_solver == "amg":
            raise ValueError("The eigen_solver was set to 'amg', but pyamg is "
                             "not available.")

    if eigen_solver is None:
        eigen_solver = 'arpack'
    elif not eigen_solver in ('arpack', 'lobpcg', 'amg'):
        raise ValueError("Unknown value for eigen_solver: '%s'."
                         "Should be 'amg', 'arpack', or 'lobpcg'"
                         % eigen_solver)

    random_state = check_random_state(random_state)

    n_nodes = adjacency.shape[0]
    # Check that the matrices given is symmetric
    if ((not sparse.isspmatrix(adjacency) and
         not np.all((adjacency - adjacency.T) < 1e-10)) or
        (sparse.isspmatrix(adjacency) and
         not np.all((adjacency - adjacency.T).data < 1e-10))):
        warnings.warn("Graph adjacency matrix should be symmetric. "
                      "Converted to be symmetric by average with its "
                      "transpose.")
    adjacency = .5 * (adjacency + adjacency.T)

    if not _graph_is_connected(adjacency):
        warnings.warn("Graph is not fully connected, spectral embedding"
                      " may not work as expected.")

    laplacian, dd = graph_laplacian(adjacency,
                                    normed=norm_laplacian, return_diag=True)

    if (eigen_solver == 'arpack'
        or eigen_solver != 'lobpcg' and
            (not sparse.isspmatrix(laplacian)
             or n_nodes < 5 * n_components)):
        # lobpcg used with eigen_solver='amg' has bugs for low number of nodes
        # for details see the source code in scipy:
        # https://github.com/scipy/scipy/blob/v0.11.0/scipy/sparse/linalg/eigen
        # /lobpcg/lobpcg.py#L237
        # or matlab:
        # http://www.mathworks.com/matlabcentral/fileexchange/48-lobpcg-m
        laplacian = _set_diag(laplacian, 1)

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
                                           sigma=1.0, which='LM',
                                           tol=eigen_tol)
            vectors = diffusion_map.T[n_components::-1]
            lambdas = lambdas[n_components::-1]
        except RuntimeError:
            # When submatrices are exactly singular, an LU decomposition
            # in arpack fails. We fallback to lobpcg
            eigen_solver = "lobpcg"

    if eigen_solver == 'amg':
        # Use AMG to get a preconditioner and speed up the eigenvalue
        # problem.
        if not sparse.issparse(laplacian):
            warnings.warn("AMG works better for sparse matrices")
        laplacian = laplacian.astype(np.float)  # lobpcg needs native floats
        laplacian = _set_diag(laplacian, 1)
        ml = smoothed_aggregation_solver(check_array(laplacian, 'csr'))
        M = ml.aspreconditioner()
        X = random_state.rand(laplacian.shape[0], n_components + 1)
        X[:, 0] = dd.ravel()
        lambdas, diffusion_map = lobpcg(laplacian, X, M=M, tol=1.e-12,
                                        largest=False)
        vectors = diffusion_map.T
        if vectors.shape[0] == 1:
            raise ValueError

    elif eigen_solver == "lobpcg":
        laplacian = laplacian.astype(np.float)  # lobpcg needs native floats
        if n_nodes < 5 * n_components + 1:
            # see note above under arpack why lobpcg has problems with small
            # number of nodes
            # lobpcg will fallback to eigh, so we short circuit it
            if sparse.isspmatrix(laplacian):
                laplacian = laplacian.todense()
            lambdas, diffusion_map = eigh(laplacian)
            vectors = diffusion_map.T[:n_components]
        else:
            # lobpcg needs native floats
            laplacian = laplacian.astype(np.float)
            laplacian = _set_diag(laplacian, 1)
            # We increase the number of eigenvectors requested, as lobpcg
            # doesn't behave well in low dimension
            X = random_state.rand(laplacian.shape[0], n_components + 1)
            X[:, 0] = dd.ravel()
            lambdas, diffusion_map = lobpcg(laplacian, X, tol=1e-15,
                                            largest=False, maxiter=2000)
            vectors = diffusion_map.T[:n_components]
            if vectors.shape[0] == 1:
                raise ValueError
    return lambdas[:n_components], vectors, dd


def spectral_embedding(adjacency, n_components=8, eigen_solver=None,
                       random_state=None, eigen_tol=0.0,
                       norm_laplacian=True, drop_first=True,
                       mode=None, diffusion_time=None):
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
    ----------
    adjacency : array-like or sparse matrix, shape: (n_samples, n_samples)
        The adjacency matrix of the graph to embed.

    n_components : integer, optional
        The dimension of the projection subspace.

    eigen_solver : {None, 'arpack', 'lobpcg', or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities.

    random_state : int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization of the
        lobpcg eigen vectors decomposition when eigen_solver == 'amg'.
        By default, arpack is used.

    eigen_tol : float, optional, default=0.0
        Stopping criterion for eigendecomposition of the Laplacian matrix
        when using arpack eigen_solver.

    drop_first : bool, optional, default=True
        Whether to drop the first eigenvector. For spectral embedding, this
        should be True as the first eigenvector should be constant vector for
        connected graph, but for spectral clustering, this should be kept as
        False to retain the first eigenvector.

    diffusion_time: float, optional, default=None
        Determines the scaling of the eigenvalues of the Laplacian

    Returns
    -------
    embedding : array, shape=(n_samples, n_components)
        The reduced samples.

    Notes
    -----
    Spectral embedding is most useful when the graph has one connected
    component. If there graph has many components, the first few eigenvectors
    will simply uncover the connected components of the graph.

    References
    ----------
    * http://en.wikipedia.org/wiki/LOBPCG

    * Toward the Optimal Preconditioned Eigensolver: Locally Optimal
      Block Preconditioned Conjugate Gradient Method
      Andrew V. Knyazev
      http://dx.doi.org/10.1137%2FS1064827500366124
    """
    # Whether to drop the first eigenvector
    if drop_first:
        n_components += 1
    _, vectors, dd = _solve_eigenvalue_problem(adjacency=adjacency,
                                               n_components=n_components,
                                               eigen_solver=eigen_solver,
                                               random_state=random_state,
                                               eigen_tol=eigen_tol,
                                               norm_laplacian=norm_laplacian,
                                               mode=mode)
    embedding = vectors * dd
    if drop_first:
        return embedding[1:].T
    else:
        return embedding.T


class SpectralEmbedding(BaseEstimator):
    """Spectral embedding for non-linear dimensionality reduction.

    Forms an affinity matrix given by the specified function and
    applies spectral decomposition to the corresponding graph laplacian.
    The resulting transformation is given by the value of the
    eigenvectors for each data point.

    Parameters
    -----------
    n_components : integer, default: 2
        The dimension of the projected subspace.

    eigen_solver : {None, 'arpack', 'lobpcg', or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities.

    random_state : int seed, RandomState instance, or None, default : None
        A pseudo random number generator used for the initialization of the
        lobpcg eigen vectors decomposition when eigen_solver == 'amg'.

    affinity : string or callable, default : "nearest_neighbors"
        How to construct the affinity matrix.
         - 'nearest_neighbors' : construct affinity matrix by knn graph
         - 'rbf' : construct affinity matrix by rbf kernel
         - 'precomputed' : interpret X as precomputed affinity matrix
         - callable : use passed in function as affinity
           the function takes in data matrix (n_samples, n_features)
           and return affinity matrix (n_samples, n_samples).

    gamma : float, optional, default : 1/n_features
        Kernel coefficient for rbf kernel.

    n_neighbors : int, default : max(n_samples/10 , 1)
        Number of nearest neighbors for nearest_neighbors graph building.

    Attributes
    ----------

    embedding_ : array, shape = (n_samples, n_components)
        Spectral embedding of the training matrix.

    affinity_matrix_ : array, shape = (n_samples, n_samples)
        Affinity_matrix constructed from samples or precomputed.

    References
    ----------

    - A Tutorial on Spectral Clustering, 2007
      Ulrike von Luxburg
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.165.9323

    - On Spectral Clustering: Analysis and an algorithm, 2011
      Andrew Y. Ng, Michael I. Jordan, Yair Weiss
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.19.8100

    - Normalized cuts and image segmentation, 2000
      Jianbo Shi, Jitendra Malik
      http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.160.2324
    """

    def __init__(self, n_components=2, affinity="nearest_neighbors",
                 gamma=None, random_state=None, eigen_solver=None,
                 n_neighbors=None):
        self.n_components = n_components
        self.affinity = affinity
        self.gamma = gamma
        self.random_state = random_state
        self.eigen_solver = eigen_solver
        self.n_neighbors = n_neighbors

    @property
    def _pairwise(self):
        return self.affinity == "precomputed"

    def _get_affinity_matrix(self, X, Y=None):
        """Caclulate the affinity matrix from data
        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

            If affinity is "precomputed"
            X : array-like, shape (n_samples, n_samples),
            Interpret X as precomputed adjacency graph computed from
            samples.

        Returns
        -------
        affinity_matrix, shape (n_samples, n_samples)
        """
        if self.affinity == 'precomputed':
            self.affinity_matrix_ = X
            return self.affinity_matrix_
        if self.affinity == 'nearest_neighbors':
            if sparse.issparse(X):
                warnings.warn("Nearest neighbors affinity currently does "
                              "not support sparse input, falling back to "
                              "rbf affinity")
                self.affinity = "rbf"
            else:
                self.n_neighbors_ = (self.n_neighbors
                                     if self.n_neighbors is not None
                                     else max(int(X.shape[0] / 10), 1))
                self.affinity_matrix_ = kneighbors_graph(X, self.n_neighbors_)
                # currently only symmetric affinity_matrix supported
                self.affinity_matrix_ = 0.5 * (self.affinity_matrix_ +
                                               self.affinity_matrix_.T)
                return self.affinity_matrix_
        if self.affinity == 'rbf':
            self.gamma_ = (self.gamma
                           if self.gamma is not None else 1.0 / X.shape[1])
            self.affinity_matrix_ = rbf_kernel(X, gamma=self.gamma_)
            return self.affinity_matrix_
        self.affinity_matrix_ = self.affinity(X)
        return self.affinity_matrix_

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

            If affinity is "precomputed"
            X : array-like, shape (n_samples, n_samples),
            Interpret X as precomputed adjacency graph computed from
            samples.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        random_state = check_random_state(self.random_state)
        if isinstance(self.affinity, six.string_types):
            if self.affinity not in set(("nearest_neighbors", "rbf",
                                         "precomputed")):
                raise ValueError(("%s is not a valid affinity. Expected "
                                  "'precomputed', 'rbf', 'nearest_neighbors' "
                                  "or a callable.") % self.affinity)
        elif not callable(self.affinity):
            raise ValueError(("'affinity' is expected to be an an affinity "
                              "name or a callable. Got: %s") % self.affinity)

        affinity_matrix = self._get_affinity_matrix(X)
        self.embedding_ = spectral_embedding(affinity_matrix,
                                             n_components=self.n_components,
                                             eigen_solver=self.eigen_solver,
                                             random_state=random_state)
        return self

    def fit_transform(self, X, y=None):
        """Fit the model from data in X and transform X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

            If affinity is "precomputed"
            X : array-like, shape (n_samples, n_samples),
            Interpret X as precomputed adjacency graph computed from
            samples.

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """
        self.fit(X)
        return self.embedding_


def diffusion_embedding(adjacency, n_components=8, diffusion_time=None,
                        eigen_solver=None, random_state=None, eigen_tol=0.0,
                        norm_laplacian=True):
    """Project the sample on the first eigenvectors of the graph Laplacian

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
    ----------
    adjacency : array-like or sparse matrix, shape: (n_samples, n_samples)
        The adjacency matrix of the graph to embed.

    n_components : integer, optional
        The dimension of the projection subspace.

    diffusion_time: float, optional, default=None
        Determines the scaling of the eigenvalues of the Laplacian

    eigen_solver : {None, 'arpack', 'lobpcg', or 'amg'}
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities.

    random_state : int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization of the
        lobpcg eigen vectors decomposition when eigen_solver == 'amg'.
        By default, arpack is used.

    eigen_tol : float, optional, default=0.0
        Stopping criterion for eigendecomposition of the Laplacian matrix
        when using arpack eigen_solver.

    Returns
    -------
    embedding : array, shape=(n_samples, n_components)
        The reduced samples.

    Notes
    -----
    Diffusion embedding is most useful when the graph has one connected
    component. If there graph has many components, the first few eigenvectors
    will simply uncover the connected components of the graph.

    References
    ----------

    - Lafon, Stephane, and Ann B. Lee. "Diffusion maps and coarse-graining: A
      unified framework for dimensionality reduction, graph partitioning, and
      data set parameterization." Pattern Analysis and Machine Intelligence,
      IEEE Transactions on 28.9 (2006): 1393-1403.
    - Coifman, Ronald R., and Stephane Lafon. Diffusion maps. Applied and
      Computational Harmonic Analysis 21.1 (2006): 5-30.

    """

    if not _graph_is_connected(adjacency):
        warnings.warn("Graph is not fully connected, spectral embedding"
                      " may not work as expected.")
    ndim = adjacency.shape[0]
    v = np.sqrt(np.sum(adjacency, axis=1))
    A = adjacency/(v[:, None] * v[None, :])
    A = np.squeeze(A * [A > 1e-5])
    print A.shape
    print n_components
    if n_components is not None:
        lambdas, vectors = eigsh(A, k=n_components)
    else:
        lambdas, vectors = eigsh(A, k=max(2, int(np.sqrt(ndim/2))))
    lambdas = lambdas[::-1]
    vectors = vectors[:, ::-1]

    psi = vectors/vectors[:, 0][:, None]
    if diffusion_time <= 0:
        lambdas = lambdas[1:] / (1 - lambdas[1:])
    else:
        lambdas = lambdas[1:]**diffusion_time
    if n_components is None:
        lambda_ratio = lambdas/lambdas[0]
        n_components = np.amin(np.nonzero(lambda_ratio < .05)[0])
        n_components = min(n_components, ndim)
    embedding = psi[:, 1:(n_components + 1)] * lambdas[:n_components][None, :]
    return embedding


class DiffusionEmbedding(SpectralEmbedding):
    """Diffusion embedding for nonlinear dimensionality reduction

    Diffusion embedding adds an extra parameter `diffusion_time` to spectral
    embedding.

    Parameters
    -----------
    diffusion_time : float
        Determines the scaling of the eigenvalues of the Laplacian

    For all other parameters see `SpectralEmbedding`

    References
    ----------

    - Lafon, Stephane, and Ann B. Lee. "Diffusion maps and coarse-graining: A
      unified framework for dimensionality reduction, graph partitioning, and
      data set parameterization." Pattern Analysis and Machine Intelligence,
      IEEE Transactions on 28.9 (2006): 1393-1403.
    - Coifman, Ronald R., and Stephane Lafon. Diffusion maps. Applied and
      Computational Harmonic Analysis 21.1 (2006): 5-30.

    """

    def __init__(self, diffusion_time=0., n_components='auto',
                 affinity="nearest_neighbors", gamma=None, random_state=None,
                 eigen_solver=None, n_neighbors=None):
        super(DiffusionEmbedding, self).__init__(n_components=n_components,
                                                 affinity=affinity,
                                                 gamma=gamma,
                                                 random_state=random_state,
                                                 eigen_solver=eigen_solver,
                                                 n_neighbors=n_neighbors)
        self.diffusion_time = diffusion_time

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

            If affinity is "precomputed"
            X : array-like, shape (n_samples, n_samples),
            Interpret X as precomputed adjacency graph computed from
            samples.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        random_state = check_random_state(self.random_state)
        if isinstance(self.affinity, basestring):
            if self.affinity not in set(("nearest_neighbors", "rbf",
                                         "precomputed", "markov")):
                raise ValueError(("%s is not a valid affinity. Expected  'markov', "
                                  "'precomputed', 'rbf', 'nearest_neighbors' "
                                  "or a callable.") % self.affinity)
        elif not hasattr(self.affinity, "__call__"):
            raise ValueError(("'affinity' is expected to be an an affinity "
                              "name or a callable. Got: %s") % self.affinity)

        if self.affinity == 'markov':
            from ..metrics import pairwise_distances
            D = pairwise_distances(X)
            k = int(max(2, np.round(D.shape[0] * 0.01)))
            eps = 2 * np.median(np.sort(D, axis=0)[k+1, :])**2
            affinity_matrix = np.exp(-(D * D) / eps)
        else:
            affinity_matrix = self._get_affinity_matrix(X)
        self.embedding_ = diffusion_embedding(affinity_matrix,
                                              n_components=self.n_components,
                                              eigen_solver=self.eigen_solver,
                                              random_state=random_state,
                                              diffusion_time=self.diffusion_time)
        return self
