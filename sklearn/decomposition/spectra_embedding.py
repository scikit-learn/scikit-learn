"""Kernel Principal Components Analysis"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Wei LI <kuantkid@gmail.com>
# License: BSD Style.

import warnings
import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_random_state
from ..utils.graph import graph_laplacian
from ..metrics.pairwise import rbf_kernel
from ..neighbors import kneighbors_graph
from ..metrics.pairwise import pairwise_kernels

from IPython.core.debugger import Tracer; debug_here = Tracer()

class SpectralEmbedding(BaseEstimator, TransformerMixin):
    """Spectral Embedding

    Non-linear dimensionality reduction through Spectral Embedding

    Project the sample on the first eigen vectors of the graph Laplacian

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

    kernel: "rbf" | "sigmoid" | "precomputed"
        Kernel.
        Default: "rbf"

    degree : int, optional
        Degree for poly, rbf and sigmoid kernels.
        Default: 3.

    gamma : float, optional
        Kernel coefficient for rbf and poly kernels.
        Default: 1/n_features.

    coef0 : float, optional
        Independent term in poly and sigmoid kernels.
    

    Attributes
    ----------

    `embedding_`
        Spectral embedding of the training matrix

    References
    ----------
    Kernel PCA was intoduced in:
        Bernhard Schoelkopf, Alexander J. Smola,
        and Klaus-Robert Mueller. 1999. Kernel principal
        component analysis. In Advances in kernel methods,
        MIT Press, Cambridge, MA, USA 327-352.
    """

    def __init__(self, n_components=None, kernel="rbf", gamma=0, degree=3,
                 coef0=1, fit_inverse_transform=False,
                 random_state=None, mode = None):
        if kernel not in {'precomputed', 'rbf', 'knn'}:
            raise ValueError(
                "Only precomputed, rbf, knn kernel supported.")
        elif fit_inverse_transform and kernel == 'precomputed':
            raise ValueError(
                "Cannot fit_inverse_transform with a precomputed kernel.")
        self.n_components = n_components
        self.kernel = kernel.lower()
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.fit_inverse_transform = fit_inverse_transform
        self.random_state = check_random_state(random_state)
        self.mode = None;

    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def _get_kernel(self, X, Y=None):
        params = {"gamma": self.gamma,
                  "degree": self.degree,
                  "coef0": self.coef0}
        if kernel == 'knn':
            return self.nn_fit.kneighbors_graph(self.nn_fit._fit_X,
            self.n_neighbors, mode='connectivity')
            
        try:
            return pairwise_kernels(X, Y, metric=self.kernel,
                                    filter_params=True, **params)
        except AttributeError:
            raise ValueError("%s is not a valid kernel. Valid kernels are: "
                             "rbf, poly, sigmoid, linear and precomputed."
                             % self.kernel)

    def _fit_transform(self, affinity):
        """ Fit's using kernel K"""
        # inverse tranform not supported
        raise NotImplementedError(
            "inverse transformation is not supported")

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """

        # get the adjacency matrix
        adjacency = self._get_kernel(X)
        # get the embedding
        self.embedding_ = self._spectra_embedding(adjacency)
        return self

    def fit_transform(self, X, y=None, **params):
        """Fit the model from data in X and transform X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """

        self.fit(X, **params)
        return self.embedding_

    def transform(self, X):
        """Transform X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """

        # out of sample not implemented
        raise NotImplementedError(
            "out of sample extension is currently unavailable")

    def _spectra_embedding(self, adjacency):
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
        except ImportError:
            if self.mode == "amg":
                raise ValueError("The mode was set to 'amg', but pyamg is "
                                 "not available.")
        n_nodes = adjacency.shape[0]
        # XXX: Should we check that the matrices given is symmetric
        if self.mode is None:
            self.mode = 'arpack'
        laplacian, dd = graph_laplacian(adjacency,
                                        normed=True, return_diag=True)
        debug_here()
        if (self.mode == 'arpack'
            or not sparse.isspmatrix(laplacian)
            or n_nodes < 5 * self.n_components):
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
            lambdas, diffusion_map = eigsh(-laplacian, k=self.n_components,
                                           sigma=1.0, which='LM')
            embedding = diffusion_map.T[::-1] * dd
        elif self.mode == 'amg':
            # Use AMG to get a preconditioner and speed up the eigenvalue
            # problem.
            laplacian = laplacian.astype(np.float)  # lobpcg needs native float
            ml = smoothed_aggregation_solver(laplacian.tocsr())
            X = self.random_state.rand(laplacian.shape[0],
                                       self.n_components)
            X[:, 0] = 1. / dd.ravel()
            M = ml.aspreconditioner()
            lambdas, diffusion_map = lobpcg(laplacian, X, M=M, tol=1.e-12,
                                            largest=False)
            embedding = diffusion_map.T * dd
            if embedding.shape[0] == 1:
                raise ValueError
        else:
            raise ValueError("Unknown value for mode: '%s'."
                             "Should be 'amg' or 'arpack'" % self.mode)
        return embedding
