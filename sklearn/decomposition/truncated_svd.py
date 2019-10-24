"""Truncated SVD for sparse matrices, aka latent semantic analysis (LSA).
"""

# Author: Lars Buitinck
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Michael Becker <mike@beckerfuffle.com>
#         Andrew Knyazev <andrew.knyazev@ucdenver.edu>
# License: 3-clause BSD.

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import svds

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array, check_random_state
from ..utils.extmath import randomized_svd
from ..utils.extmath import safe_sparse_dot, svd_flip
from ..utils.sparsefuncs import mean_variance_axis

__all__ = ["TruncatedSVD"]


class TruncatedSVD(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using truncated SVD (aka LSA).

    This transformer performs linear dimensionality reduction by means of
    truncated singular value decomposition (SVD). Contrary to PCA, this
    estimator does not center the data before computing the singular value
    decomposition. This means it can work with scipy.sparse matrices
    efficiently.

    In particular, truncated SVD works on term count/tf-idf matrices as
    returned by the vectorizers in sklearn.feature_extraction.text. In that
    context, it is known as latent semantic analysis (LSA).

    This estimator supports 3 algorithms: ARPACK (Lanczos), randomized
    (block power), and LOBPCG (block conjugate gradient type) iterations, used
    within a "naive" approach that reduces SVD to the eigenproblem for
    the normal matrix (X * X.T) or (X.T * X), whichever has a smaller size.

    Read more in the :ref:`User Guide <LSA>`.

    Parameters
    ----------
    n_components : int, default=2
        Desired dimensionality of output data.
        Must be strictly less than the number of features.
        The default value is useful for visualisation. For LSA, a value of
        100 is recommended.

    algorithm : string {'randomized', 'lobpcg', 'arpack'}, default="randomized"
        SVD solver to use.
        randomized :
            run randomized SVD due to Halko (2009).
        lobpcg :
            run Locally Optimal Block Preconditioned Conjugate Gradient [2]
            for a normal matrix X'*X or X*X', whichever of the two is of
            the smallest size. See :func:`scipy.sparse.linalg.lobpcg`.
        arpack :
            run SVD truncated to n_components calling ARPACK solver via
            :func:`scipy.sparse.linalg.svds`. It requires strictly
            0 < n_components < min(X.shape).

    n_iter : int, default=5
        Number of iterations for algorithm 'randomized' and 'lobpcg'. Not used
        by ARPACK.
        The default is larger than the default in
        `~sklearn.utils.extmath.randomized_svd` to handle
        sparse matrices that may have large slowly decaying spectrum.

    random_state : int, RandomState instance or None, optional, default = None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    tol : None or float, default=None
        Tolerance for singular values computed by the 'randomized' and 'lobpcg'
        SVD solver. If None, then:
            * `tol = 2 * eps` for svd_solver = 'randomized'.
              Refer to :func:`scipy.sparse.linalg.svds`.
            * `tol = n * sqrt(eps)` where `n = min(n_samples, n_features)`.
              Refer to :func:`scipy.sparse.linalg.lobpcg`.

    Attributes
    ----------
    components_ : array, shape (n_components, n_features)

    explained_variance_ : array, shape (n_components,)
        The variance of the training samples transformed by a projection to
        each component.

    explained_variance_ratio_ : array, shape (n_components,)
        Percentage of variance explained by each of the selected components.

    singular_values_ : array, shape (n_components,)
        The singular values corresponding to each of the selected components.
        The singular values are equal to the 2-norms of the ``n_components``
        variables in the lower-dimensional space.

    References
    ----------
    .. [1] Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp.
            "Finding structure with randomness: Probabilistic algorithms for
            constructing approximate matrix decompositions."
            SIAM review 53, no. 2 (2011): 217-288.

    .. [2] Knyazev, Andrew V.
            "Toward the optimal preconditioned eigensolver:
            Locally optimal block preconditioned conjugate gradient method."
            SIAM journal on scientific computing 23, no. 2 (2001): 517-541.
            https://doi.org/10.1137%2FS1064827500366124

    Examples
    --------
    >>> from sklearn.decomposition import TruncatedSVD
    >>> from sklearn.random_projection import sparse_random_matrix
    >>> X = sparse_random_matrix(100, 100, density=0.01, random_state=42)
    >>> svd = TruncatedSVD(n_components=5, n_iter=7, random_state=42)
    >>> svd.fit(X)
    TruncatedSVD(n_components=5, n_iter=7, random_state=42)
    >>> print(svd.explained_variance_ratio_)
    [0.0606... 0.0584... 0.0497... 0.0434... 0.0372...]
    >>> print(svd.explained_variance_ratio_.sum())
    0.249...
    >>> print(svd.singular_values_)
    [2.5841... 2.5245... 2.3201... 2.1753... 2.0443...]

    See also
    --------
    PCA

    Notes
    -----
    SVD suffers from a problem called "sign indeterminacy", which means the
    sign of the ``components_`` and the output from transform depend on the
    algorithm and random state. To work around this, fit instances of this
    class to data once, then keep the instance around to do transformations.

    """

    def __init__(self, n_components=2, algorithm="randomized", n_iter=5,
                 random_state=None, tol=None):
        self.algorithm = algorithm
        self.n_components = n_components
        self.n_iter = n_iter
        self.random_state = random_state
        self.tol = tol

    def fit(self, X, y=None):
        """Fit LSI model on training data X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data.

        y : Ignored

        Returns
        -------
        self : object
            Returns the transformer object.
        """
        self.fit_transform(X)
        return self

    def fit_transform(self, X, y=None):
        """Fit LSI model to X and perform dimensionality reduction on X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data.

        y : Ignored

        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Reduced version of X. This will always be a dense array.
        """
        X = check_array(X, accept_sparse=['csr', 'csc'],
                        ensure_min_features=2)
        random_state = check_random_state(self.random_state)

        tol = (0. if self.tol is None and self.algorithm == 'arpack'
               else self.tol)

        if self.algorithm == "arpack":
            U, Sigma, VT = svds(X, k=self.n_components, tol=tol)
            # svds doesn't abide by scipy.linalg.svd/randomized_svd
            # conventions, so reverse its outputs.
            Sigma = Sigma[::-1]
            U, VT = svd_flip(U[:, ::-1], VT[::-1])
        elif self.algorithm == "lobpcg":
            k = self.n_components
            n_features = X.shape[1]
            if k >= n_features:
                raise ValueError("n_components must be < n_features;"
                                 " got %d >= %d" % (k, n_features))
            U, Sigma, VT = randomized_svd(
                X, self.n_components, n_iter=self.n_iter,
                random_state=random_state,
                preconditioner='lobpcg', tol=tol
            )
        elif self.algorithm == "randomized":
            k = self.n_components
            n_features = X.shape[1]
            if k >= n_features:
                raise ValueError("n_components must be < n_features;"
                                 " got %d >= %d" % (k, n_features))
            U, Sigma, VT = randomized_svd(X, self.n_components,
                                          n_iter=self.n_iter,
                                          random_state=random_state)
        else:
            raise ValueError("unknown algorithm %r" % self.algorithm)

        self.components_ = VT

        # Calculate explained variance & explained variance ratio
        X_transformed = U * Sigma
        self.explained_variance_ = exp_var = np.var(X_transformed, axis=0)
        if sp.issparse(X):
            _, full_var = mean_variance_axis(X, axis=0)
            full_var = full_var.sum()
        else:
            full_var = np.var(X, axis=0).sum()
        self.explained_variance_ratio_ = exp_var / full_var
        self.singular_values_ = Sigma  # Store the singular values.

        return X_transformed

    def transform(self, X):
        """Perform dimensionality reduction on X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            New data.

        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Reduced version of X. This will always be a dense array.
        """
        X = check_array(X, accept_sparse='csr')
        return safe_sparse_dot(X, self.components_.T)

    def inverse_transform(self, X):
        """Transform X back to its original space.

        Returns an array X_original whose transform would be X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_components)
            New data.

        Returns
        -------
        X_original : array, shape (n_samples, n_features)
            Note that this is always a dense array.
        """
        X = check_array(X)
        return np.dot(X, self.components_)
