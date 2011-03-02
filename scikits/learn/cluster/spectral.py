""" Algorithms for  spectral clustering.
"""

# Author: Gael Varoquaux gael.varoquaux@normalesup.org
# License: BSD
import warnings

import numpy as np


from ..base import BaseEstimator
from ..utils.graph import graph_laplacian
from .k_means_ import k_means


def spectral_embedding(adjacency, k=8, mode=None):
    """ Spectral embedding: project the sample on the k first
        eigen vectors of the normalized graph Laplacian.

        Parameters
        -----------
        adjacency: array-like or sparse matrix, shape: (p, p)
            The adjacency matrix of the graph to embed.
        k: integer, optional
            The dimension of the projection subspace.
        mode: {None, 'arpack' or 'amg'}
            The eigenvalue decomposition strategy to use. AMG (Algebraic
            MultiGrid) is much faster, but requires pyamg to be
            installed.

        Returns
        --------
        embedding: array, shape: (p, k)
            The reduced samples

        Notes
        ------
        The graph should contain only one connect component,
        elsewhere the results make little sens.
    """

    from scipy import sparse
    from ..utils.fixes import arpack_eigsh
    from scipy.sparse.linalg import lobpcg
    try:
        from pyamg import smoothed_aggregation_solver
        amg_loaded = True
    except ImportError:
        amg_loaded = False

    n_nodes = adjacency.shape[0]
    # XXX: Should we check that the matrices given is symmetric
    if not amg_loaded:
        warnings.warn('pyamg not available, using scipy.sparse')
    if mode is None:
        mode = ('amg' if amg_loaded else 'arpack')
    laplacian, dd = graph_laplacian(adjacency,
                                    normed=True, return_diag=True)
    if (mode == 'arpack'
        or not sparse.isspmatrix(laplacian)
        or n_nodes < 5*k # This is the threshold under which lobpcg has bugs
       ):
        # We need to put the diagonal at zero
        if not sparse.isspmatrix(laplacian):
            laplacian[::n_nodes+1] = 0
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
        lambdas, diffusion_map = arpack_eigsh(-laplacian, k=k, which='LA')
        embedding = diffusion_map.T[::-1]*dd
    elif mode == 'amg':
        # Use AMG to get a preconditionner and speed up the eigenvalue
        # problem.
        laplacian = laplacian.astype(np.float) # lobpcg needs the native float
        ml = smoothed_aggregation_solver(laplacian.tocsr())
        X = np.random.rand(laplacian.shape[0], k)
        X[:, 0] = 1. / dd.ravel()
        M = ml.aspreconditioner()
        lambdas, diffusion_map = lobpcg(laplacian, X, M=M, tol=1.e-12,
                                        largest=False)
        embedding = diffusion_map.T * dd
        if embedding.shape[0] == 1: raise ValueError
    else:
        raise ValueError("Unknown value for mode: '%s'." % mode)
    return embedding


def spectral_clustering(adjacency, k=8, mode=None):
    """ Spectral clustering: apply k-means to a projection of the
        graph laplacian, finds normalized graph cuts.

        Parameters
        -----------
        adjacency: array-like or sparse matrix, shape: (p, p)
            The adjacency matrix of the graph to embed.
        k: integer, optional
            The dimension of the projection subspace.
        mode: {None, 'arpack' or 'amg'}
            The eigenvalue decomposition strategy to use. AMG (Algebraic
            MultiGrid) is much faster, but requires pyamg to be
            installed.

        Returns
        --------
        labels: array of integers, shape: p
            The labels of the clusters.
        centers: array of integers, shape: k
            The indices of the cluster centers

        Notes
        ------
        The graph should contain only one connect component,
        elsewhere the results make little sens.

        This algorithm solves the normalized cut for k=2: it is a
        normalized spectral clustering.
    """
    maps = spectral_embedding(adjacency, k=k, mode=mode)
    maps = maps[1:]
    _, labels, _ = k_means(maps.T, k)
    return labels


################################################################################
class SpectralClustering(BaseEstimator):
    """ Spectral clustering: apply k-means to a projection of the
        graph laplacian, finds normalized graph cuts.

        Parameters
        -----------
        k: integer, optional
            The dimension of the projection subspace.
        mode: {None, 'arpack' or 'amg'}
            The eigenvalue decomposition strategy to use. AMG (Algebraic
            MultiGrid) is much faster, but requires pyamg to be
            installed.

        Methods
        -------

        fit(X):
            Compute spectral clustering

        Attributes
        ----------

        labels_:
            Labels of each point

    """


    def __init__(self, k=8, mode=None):
        self.k = k
        self.mode = mode


    def fit(self, X, **params):
        """ Compute the spectral clustering from the adjacency matrix of
            the graph.

            Parameters
            -----------
            X: array-like or sparse matrix, shape: (p, p)
                The adjacency matrix of the graph to embed.
                X is an adjacency matrix of a similarity graph: its
                entries must be positive or zero. Zero means that
                elements have nothing in common, whereas high values mean
                that elements are strongly similar.

            Notes
            ------
            If you have an affinity matrix, such as a distance matrix,
            for which 0 means identical elements, and high values means
            very dissimilar elements, it can be transformed in a
            similarity matrix that is well suited for the algorithm by
            applying the gaussian (heat) kernel::

                np.exp(- X**2/2. * delta**2)

            If the pyamg package is installed, it is used. This
            greatly speeds up computation.
        """
        self._set_params(**params)
        self.labels_ = spectral_clustering(X,
                                k=self.k, mode=self.mode)
        return self

