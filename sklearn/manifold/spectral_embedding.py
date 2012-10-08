"""Spectral Embedding"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Wei LI <kuantkid@gmail.com>
# License: BSD Style.

import warnings
import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_random_state
from ..utils.graph import graph_laplacian
from ..utils.arpack import eigsh
from ..metrics.pairwise import rbf_kernel
from ..neighbors import kneighbors_graph
from ..metrics.pairwise import pairwise_kernels
from ..metrics.pairwise import rbf_kernel
from scipy import sparse
from scipy.sparse.linalg import lobpcg


def spectral_embedding(adjacency, n_components=2,
                       eigen_solver=None, random_state=None,
                       norm_laplacian=True):
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
    adjacency : array-like or sparse matrix, shape: (n_samples, n_samples)
        The adjacency matrix of the graph to embed.

    n_components : integer, default: 2
        The dimension of the projected subspace.

    norm_laplacian : bool, default: False
        Use normalized graph laplacian (True) or un-normalized graph
        laplacian (False)

    Returns
    --------
    embedding : array, shape: (n_samples, n_components)
        The reduced samples.

    Notes
    ------
    The graph should contain only one connected component, elsewhere the
    results make little sense.
    """

    try:
        from pyamg import smoothed_aggregation_solver
    except ImportError:
        if eigen_solver == "amg":
            raise ValueError(
                "The eigen_solver was set to 'amg', but pyamg is "
                "not available.")
    n_nodes = adjacency.shape[0]
    # XXX: Should we check that the matrices given is symmetric
    if eigen_solver is None:
        eigen_solver = 'arpack'
    laplacian, dd = graph_laplacian(adjacency,
                                    normed=True, return_diag=True)
    if (eigen_solver == 'arpack'
        or not sparse.isspmatrix(laplacian)
        or n_nodes < 5 * n_components):
        # lobpcg used with mode='amg' has bugs for low number of nodes

        # We need to put the diagonal at zero
        if not sparse.isspmatrix(laplacian):
            laplacian.flat[::n_nodes + 1] = 0
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

        # Here we'll use shift-invert mode for fast eigenvalues (see
        # http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html
        # for a short explanation of what this means) Because the
        # normalized Laplacian has eigenvalues between 0 and 2, I - L has
        # eigenvalues between -1 and 1.  ARPACK is most efficient when
        # finding eigenvalues of largest magnitude (keyword which='LM') and
        # when these eigenvalues are very large compared to the rest. For
        # very large, very sparse graphs, I - L can have many, many
        # eigenvalues very near 1.0.  This leads to slow convergence.  So
        # instead, we'll use ARPACK's shift-invert mode, asking for the
        # eigenvalues near 1.0.  This effectively spreads-out the spectrum
        # near 1.0 and leads to much faster convergence: potentially an
        # orders-of-magnitude speedup over simply using keyword which='LA'
        # in standard mode.
        lambdas, diffusion_map = eigsh(-laplacian, k=n_components + 1,
                                        sigma=1.0, which='LM')
        embedding = diffusion_map.T[n_components - 1::-1] * dd
    elif eigen_solver == 'amg':
        # Use AMG to get a preconditioner and speed up the eigenvalue
        # problem.
        laplacian = laplacian.astype(np.float)  # lobpcg needs native float
        ml = smoothed_aggregation_solver(laplacian.tocsr())
        random_state = check_random_state(random_state)
        X = random_state.rand(laplacian.shape[0],
                                   n_components)
        X[:, 0] = 1. / dd.ravel()
        M = ml.aspreconditioner()
        lambdas, diffusion_map = lobpcg(laplacian, X, M=M, tol=1.e-12,
                                        largest=False)
        embedding = diffusion_map.T * dd
        if embedding.shape[0] == 1:
            raise ValueError
    else:
        raise ValueError("Unknown value for mode: '%s'."
                         "Should be 'amg' or 'arpack'" % eigen_solver)
    return embedding.T


class SpectralEmbedding(BaseEstimator, TransformerMixin):
    """Spectral Embedding for Non-linear Dimensionality Reduction.

    Forms an affinity matrix given by the specified function and
    applies spectral decomposition to the corresponding graph laplacian.
    The resulting transformation is given by the value of the
    eigenvectors for each data point.

    Parameters
    -----------
    n_components : integer, default: 2
        The dimension of the projected subspace.

    eigen_solver : {None, 'arpack' or 'amg'}
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

    fit_inverse_transform : bool, optional, default : False
       Whether to fit the inverse transformation.

    Attributes
    ----------

    `embedding_` : array, shape = (n_samples, n_components)
        Spectral embedding of the training matrix.

    `affinity_matrix_` : array, shape = (n_samples, n_samples)
        Affinity_matrix constructed from samples or precomputed.

    References
    ----------
    Spectral Embedding was intoduced in:
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
                 gamma=None, fit_inverse_transform=False,
                 random_state=None, eigen_solver=None, n_neighbors=None):
        self.n_components = n_components
        self.affinity = affinity
        self.gamma = gamma
        self.fit_inverse_transform = fit_inverse_transform
        self.random_state = check_random_state(random_state)
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
            if self.gamma is None:
                self.gamma = 1.0 / X.shape[1]
            if self.n_neighbors is None:
                self.n_neighbors = max(int(X.shape[0] / 10), 1)
            self.affinity_matrix_ = kneighbors_graph(X, self.n_neighbors)
            # currently only symmetric affinity_matrix supported
            self.affinity_matrix_ = 0.5 * (self.affinity_matrix_ +
                                           self.affinity_matrix_.T)
            return self.affinity_matrix_
        if self.affinity == 'rbf':
            if self.gamma is None:
                self.gamma = 1.0 / X.shape[1]
            self.affinity_matrix_ = rbf_kernel(X, gamma=self.gamma)
            return self.affinity_matrix_
        try:
            self.affinity_matrix_ = self.affinity(X)
            return self.affinity_matrix_
        except:
            raise ValueError(
                "%s is not a valid graph type. Valid kernels are: "
                "nearest_neighbors, precomputed and callable."
                % self.affinity)

    def fit(self, X, y=None):
        """Fit the model from data in X.

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
        self : object
            Returns the instance itself.
        """
        if isinstance(self.affinity, str):
            if self.affinity not in \
                {'precomputed', 'rbf', 'nearest_neighbors'}:
                raise ValueError(
                    "Only precomputed, rbf,"
                    "nearest_neighbors graph supported.")
        if self.fit_inverse_transform and self.affinity == 'precomputed':
            raise ValueError(
                "Cannot fit_inverse_transform with a precomputed kernel.")
        affinity_matrix = self._get_affinity_matrix(X)
        self.embedding_ = spectral_embedding(affinity_matrix,
                                             n_components=self.n_components,
                                             eigen_solver=self.eigen_solver,
                                             random_state=self.random_state)
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
