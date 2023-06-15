"""Graph Spectral Embedding."""

# Author: Ali Saad-Eldin <ali.saadeldin11@gmail.com>

# License: BSD 3 clause

import warnings

import numpy as np
from scipy import sparse
from scipy.linalg import svd
from scipy.sparse import diags, isspmatrix_csr
from scipy.stats import norm

from scipy.sparse.csgraph import connected_components
from ..base import BaseEstimator
from ..utils import check_symmetric
from ..utils.extmath import randomized_svd
from ..utils.validation import check_is_fitted


def _augment_diagonal(graph, weight=1):
    r"""
    Replaces the diagonal of an adjacency matrix with :math:`\frac{d}{nverts - 1}` where
    :math:`d` is the degree vector for an unweighted graph and the sum of magnitude of
    edge weights for each node for a weighted graph. For a directed graph the in/out
    :math:`d` is averaged.

    Parameters
    ----------
    graph: {ndarray, sparse matrix}
        Input graph in any of the above specified formats.

    weight: float/int
        scalar value to multiply the new diagonal vector by

    Returns
    -------
    graph : np.array
        Adjacency matrix with average degrees added to the diagonal.
    """
    dia = diags(graph.diagonal()) if isspmatrix_csr(graph) else np.diag(np.diag(graph))
    graph = graph - dia

    divisor = graph.shape[0] - 1

    in_degrees = np.squeeze(np.asarray(abs(graph).sum(axis=0)))
    out_degrees = np.squeeze(np.asarray(abs(graph).sum(axis=1)))

    degrees = (in_degrees + out_degrees) / 2
    diag = weight * degrees / divisor

    graph += diags(diag) if isspmatrix_csr(graph) else np.diag(diag)

    return graph


def _selectSVD(X, n_components=None, n_elbows=2, svd_solver="randomized"):
    r"""
    Dimensionality reduction using either full Singular Value Decomposition
    (SVD) or truncated SVD.

    Full SVD is performed using SciPy's wrapper for ARPACK,
    while truncated SVD is performed using either SciPy's wrapper for LAPACK
    or scikit-learn's implementation of randomized SVD.

    It also performs optimal dimensionality selection using Zhu & Godsie [1]_
    `svd_solver` if the number of target dimension is not specified.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        The data to perform whose dimension are to reduce.

    n_components : int or None, default=None
        Desired dimensionality of output data. If "full",
        ``n_components`` must be ``<= min(X.shape[0])``.
        Otherwise, ``n_components`` must be ``< min(X.shape[0])``.
        If `None`, then optimal dimensions will be chosen
        using ``n_elbows`` argument.

    n_elbows : int, optional, default=2
        If ``n_components`` is None, then compute the optimal embedding dimension using
        [1]_. Otherwise, ignored.

    svd_solver : {'randomized' (default), 'full', 'truncated'}, optional
        SVD solver to use:
        - 'randomized'
            Computes randomized SVD using
            :func:`sklearn.utils.extmath.randomized_svd`
        - 'full'
            Computes full SVD using :func:`scipy.linalg.svd`
            Does not support ``graph`` input of type `scipy.sparse.csr_matrix`
        - 'truncated'
            Computes truncated SVD using :func:`scipy.sparse.linalg.svds`

    Returns
    -------
    U : array-like, shape (n_samples, n_components)
        Left singular vectors corresponding to singular values.

    D : array-like, shape (n_components)
        Singular values in decreasing order, as a 1d array.

    V : array-like, shape (n_components, n_samples)
        Right singular vectors corresponding to singular values.

    References
    ----------
    .. [1] Zhu, M. and Ghodsi, A. (2006).
        Automatic dimensionality selection from the scree plot via the use of
        profile likelihood. Computational Statistics & Data Analysis, 51(2),
        pp.918-930.
    """
    # Added in order to pass check estimator, must include words "one sample"
    if X.shape[0] == 1:
        msg = "Input data has only one sample (node)"
        raise ValueError(msg)

    # Deal with svd_solvers
    if svd_solver not in ["full", "truncated", "randomized"]:
        msg = "svd_solver must be one of {full, truncated, randomized}."
        raise ValueError(msg)

    if svd_solver == "full" and isspmatrix_csr(X):
        msg = "'full' agorithm does not support scipy.sparse.csr_matrix inputs."
        raise TypeError(msg)

    if n_components is None:
        elbows, _ = _select_dimension(X, n_elbows=n_elbows, threshold=None)
        n_components = elbows[-1]
    elif not isinstance(n_components, int):
        msg = "n_components must be an integer, not {}.".format(type(n_components))
        raise ValueError(msg)

    # Check
    if (svd_solver == "full") & (n_components > min(X.shape)):
        msg = "n_components must be <= min(X.shape)."
        raise ValueError(msg)
    elif svd_solver == "full":
        U, D, V = svd(X)
        U = U[:, :n_components]
        D = D[:n_components]
        V = V[:n_components, :]

    if (svd_solver in ["truncated", "randomized"]) & (n_components >= min(X.shape)):
        msg = "n_components must be strictly < min(X.shape)."
        raise ValueError(msg)
    elif svd_solver == "truncated":
        U, D, V = sparse.linalg.svds(X, k=n_components)
        idx = np.argsort(D)[::-1]  # sort in decreasing order
        D = D[idx]
        U = U[:, idx]
        V = V[idx, :]
    elif svd_solver == "randomized":
        U, D, V = randomized_svd(X, n_components)

    return U, D, V


def _select_dimension(
    X, n_components=None, n_elbows=2, threshold=None, return_likelihoods=False
):
    """
    Generates profile likelihood from array based on Zhu and Godsie method.
    Elbows correspond to the optimal embedding dimension.

    Parameters
    ----------
    X : 1d or 2d array-like
        Input array generate profile likelihoods for. If 1d-array, it should be
        sorted in decreasing order. If 2d-array, shape should be
        (n_samples, n_features).

    n_components : int, optional, default: None.
        Number of components to embed. If None, ``n_components =
        floor(log2(min(n_samples, n_features)))``. Ignored if ``X`` is 1d-array.

    n_elbows : int, optional, default: 2.
        Number of likelihood elbows to return. Must be ``> 1``.

    threshold : float, int, optional, default: None
        If given, only consider the singular values that are ``> threshold``. Must
        be ``>= 0``.

    return_likelihoods : bool, optional, default: False
        If True, returns the all likelihoods associated with each elbow.

    Returns
    -------
    elbows : list
        Elbows indicate subsequent optimal embedding dimensions. Number of
        elbows may be less than ``n_elbows`` if there are not enough singular
        values.

    sing_vals : list
        The singular values associated with each elbow.

    likelihoods : list of array-like
        Array of likelihoods of the corresponding to each elbow. Only returned
        if ``return_likelihoods`` is True.

    References
    ----------
    .. [1] Zhu, M. and Ghodsi, A. (2006).
        Automatic dimensionality selection from the scree plot via the use of
        profile likelihood. Computational Statistics & Data Analysis, 51(2),
        pp.918-930.
    """
    # Handle input data
    if not isinstance(X, np.ndarray) and not isspmatrix_csr(X):
        msg = "X must be a numpy array or scipy.sparse.csr_matrix, not {}.".format(
            type(X)
        )
        raise ValueError(msg)
    if X.ndim > 2:
        msg = "X must be a 1d or 2d-array, not {}d array.".format(X.ndim)
        raise ValueError(msg)
    elif np.min(X.shape) <= 1:
        msg = "X must have more than 1 samples or 1 features."
        raise ValueError(msg)

    # Handle n_elbows
    if not isinstance(n_elbows, int):
        msg = "n_elbows must be an integer, not {}.".format(type(n_elbows))
        raise ValueError(msg)
    elif n_elbows < 1:
        msg = "number of elbows should be an integer > 1, not {}.".format(n_elbows)
        raise ValueError(msg)

    # Handle threshold
    if threshold is not None:
        if not isinstance(threshold, (int, float)):
            msg = "threshold must be an integer or a float, not {}.".format(
                type(threshold)
            )
            raise ValueError(msg)
        elif threshold < 0:
            msg = "threshold must be >= 0, not {}.".format(threshold)
            raise ValueError(msg)

    # Handle n_components
    if n_components is None:
        # per recommendation by Zhu & Godsie
        k = int(np.ceil(np.log2(np.min(X.shape))))
    elif not isinstance(n_components, int):
        msg = "n_components must be an integer, not {}.".format(type(n_components))
        raise ValueError(msg)
    else:
        k = n_components

    # Check to see if svd is needed
    if X.ndim == 1:
        D = np.sort(X)[::-1]
    elif X.ndim == 2:
        # Singular values in decreasing order
        D = sparse.linalg.svds(A=X, k=k, return_singular_vectors=False)
        D = np.sort(D)[::-1]
        # U, D, V = sklearn.utils.extmath.randomized_svd()

    if threshold is not None:
        D = D[D > threshold]

        if len(D) == 0:
            msg = "No values greater than threshold {}."
            raise IndexError(msg.format(threshold))

    idx = 0
    elbows = []
    values = []
    likelihoods = []
    for _ in range(n_elbows):
        arr = D[idx:]
        if arr.size <= 1:  # Cant compute likelihoods with 1 numbers
            break
        lq = _compute_likelihood(arr)
        idx += np.argmax(lq) + 1
        elbows.append(idx)
        values.append(D[idx - 1])
        likelihoods.append(lq)

    if return_likelihoods:
        return elbows, values, likelihoods
    else:
        return elbows, values


def _compute_likelihood(arr):
    """
    Computes the log likelihoods based on normal distribution given
    a 1d-array of sorted values. If the input has no variance,
    the likelihood will be nan.
    """
    n_elements = len(arr)
    likelihoods = np.zeros(n_elements)

    for idx in range(1, n_elements + 1):
        # split into two samples
        s1 = arr[:idx]
        s2 = arr[idx:]

        # deal with when input only has 2 elements
        if (s1.size == 1) & (s2.size == 1):
            likelihoods[idx - 1] = -np.inf
            continue

        # compute means
        mu1 = np.mean(s1)
        if s2.size != 0:
            mu2 = np.mean(s2)
        else:
            # Prevent numpy warning for taking mean of empty array
            mu2 = -np.inf

        # compute pooled variance
        variance = ((np.sum((s1 - mu1) ** 2) + np.sum((s2 - mu2) ** 2))) / (
            n_elements - 1 - (idx < n_elements)
        )
        std = np.sqrt(variance)

        # compute log likelihoods
        likelihoods[idx - 1] = np.sum(norm.logpdf(s1, loc=mu1, scale=std)) + np.sum(
            norm.logpdf(s2, loc=mu2, scale=std)
        )

    return likelihoods


def _to_laplacian(A, regularizer=None):
    r"""
    A function to convert graph adjacency matrix to graph Laplacian.

    Currently supports I-DAD, DAD, and R-DAD Laplacians, where D is the diagonal
    matrix of degrees of each node raised to the -1/2 power, I is the
    identity matrix, and A is the adjacency matrix.

    R-DAD is regularized Laplacian: where :math:`D_t = D + regularizer \times I`.

    Parameters
    ----------
    graph: object
        Either array-like, (n_vertices, n_vertices) numpy array,
        scipy.sparse.csr_matrix, or an object of type networkx.Graph.

    form: {'I-DAD' (default), 'DAD', 'R-DAD'}, string, optional
        - 'I-DAD'
            Computes :math:`L = I - D_i A D_i`
        - 'DAD'
            Computes :math:`L = D_o A D_i`
        - 'R-DAD'
            Computes :math:`L = D_o^r A D_i^r`
            where :math:`D_o^r = D_o + regularizer \times I` and likewise for
            :math:`D_i`

    regularizer: int, float or None, optional (default=None)
        Constant to add to the degree vector(s). If None, average node degree is added.
        If int or float, must be >= 0. Only used when ``form`` is 'R-DAD'.

    Returns
    -------
    L : numpy.ndarray
        2D (n_vertices, n_vertices) array representing graph
        Laplacian of specified form

    References
    ----------
    .. [1] Qin, Tai, and Karl Rohe. "Regularized spectral clustering
           under the degree-corrected stochastic blockmodel." In Advances
           in Neural Information Processing Systems, pp. 3120-3128. 2013
    .. [2] Rohe, Karl, Tai Qin, and Bin Yu. "Co-clustering directed graphs to discover
           asymmetries and directional communities." Proceedings of the National Academy
           of Sciences 113.45 (2016): 12679-12684.
    """

    in_degree = np.reshape(np.asarray(A.sum(axis=0)), (-1,))
    out_degree = np.reshape(np.asarray(A.sum(axis=1)), (-1,))

    # regularize laplacian with parameter
    # set to average degree
    # if regularizer == True or int:
    if isinstance(regularizer, bool) and regularizer:
        regularizer = np.mean(out_degree)
    elif not isinstance(regularizer, (int, float, bool)):
        raise TypeError(
            "Regularizer must be a int or float, not {}".format(type(regularizer))
        )
    elif regularizer < 0:
        raise ValueError("Regularizer must be greater than or equal to 0")

    if isinstance(regularizer, (int, float)):
        in_degree += regularizer
        out_degree += regularizer

    with np.errstate(divide="ignore"):
        in_root = 1 / np.sqrt(in_degree)  # this is 10x faster than ** -0.5
        out_root = 1 / np.sqrt(out_degree)

    diag = diags if isspmatrix_csr(A) else np.diag

    in_root[np.isinf(in_root)] = 0
    out_root[np.isinf(out_root)] = 0

    in_root = diag(in_root)  # just change to sparse diag for sparse support
    out_root = diag(out_root)

    L = out_root @ A @ in_root
    if _is_symmetric(A):
        return check_symmetric(L)
        # sometimes machine prec. makes this necessary
    return L


def _is_symmetric(array, tol=1e-15):
    """Check if matrix is symmetric.

    Parameters
    ----------
    array : {ndarray, sparse matrix}
        Input object to check / convert. Must be two-dimensional and square,
        otherwise a ValueError will be raised.

    tol : float, default=1e-10
        Absolute tolerance for equivalence of arrays. Default = 1E-15.

    Returns
    -------
    symmetric : bool
        Whether or not array is symmetric
    """
    if isspmatrix_csr(array):
        diff = array - array.T
        symmetric = np.all(abs(diff.data) < tol)
    else:
        symmetric = np.allclose(array, array.T, atol=tol)

    return symmetric


class GraphSpectralEmbedding(BaseEstimator):
    r"""Spectral embedding for dimensionality reduction of graph represented data.

    Computes the adjacency or laplacian spectral embedding of a graph.
    The adjacency spectral embedding (ASE) is a k-dimensional Euclidean representation
    of the graph based on its adjacency matrix, while the laplacian spectral embedding
    (LSE) is based on its Laplacian matrix. It relies on an SVD to reduce the
    dimensionality to the specified k, or if k is unspecified, can find a number of
    dimensions automatically.


    Parameters
    ----------
    n_components : int or None, default = None
        Desired dimensionality of output data. If "full",
        ``n_components`` must be ``<= min(X.shape)``. Otherwise, ``n_components`` must
        be ``< min(X.shape)``. If None, then optimal dimensions will be chosen by Zhu
        and Ghodsi method using ``n_elbows`` argument.

    n_elbows : int, optional, default: 2
        If ``n_components`` is None, then compute the optimal embedding dimension using
        Zhu and Ghodsi method. Otherwise, ignored.

    aglorithm : {'ASE' , 'LSE'}, default = 'ASE'
        - 'ASE'
            Adjacency Spectral Embedding
        - 'LSE'
            Laplacian Spectral Embedding

    svd_solver : {'randomized' (default), 'full', 'truncated'}, optional
        SVD solver to use:
        - 'randomized'
            Computes randomized svd using
            :func:`sklearn.utils.extmath.randomized_svd`
        - 'full'
            Computes full svd using :func:`scipy.linalg.svd`
            Does not support ``graph`` input of type scipy.sparse.csr_matrix
        - 'truncated'
            Computes truncated svd using :func:`scipy.sparse.linalg.svds`

    check_lcc : bool , optional (default = True)
        Whether to check if input graph is connected. May result in non-optimal
        results if the graph is unconnected. If True and input is unconnected,
        a UserWarning is thrown. Not checking for connectedness may result in
        faster computation.

    regularizer: int, float or bool, optional (default = True)
        When `algorithm`='ASE', whether to replace the main diagonal of the adjacency
        matrix with a vector corresponding to the degree (or sum of edge weights for a
        weighted network) before embedding. Empirically, this produces latent position
        estimates closer to the ground truth.

        When `algorithm`='LSE', constant to be added to the diagonal of degree matrix.
        If `True`, average node degree is added. If int or float, must be >= 0.

    concat : bool, optional (default False)
        If graph is directed, whether to concatenate left and right (out and in) latent
        positions along axis 1.

    Attributes
    ----------
    n_features_in_: int
        Number of features passed to the method.

    latent_left_ : array, shape (n_samples, n_components)
        Estimated left latent positions of the graph.

    latent_right_ : array, shape (n_samples, n_components), or None
        Only computed when the graph is directed, or adjacency matrix is assymetric.
        Estimated right latent positions of the graph. Otherwise, None.

    singular_values_ : array, shape (n_components)
        Singular values associated with the latent position matrices.

    Notes
    -----
    The singular value decomposition:
    .. math:: A = U \Sigma V^T
    is used to find an orthonormal basis for a matrix, which in our case is the
    adjacency matrix of the graph. These basis vectors (in the matrices U or V) are
    ordered according to the amount of variance they explain in the original matrix.
    By selecting a subset of these basis vectors (through our choice of dimensionality
    reduction) we can find a lower dimensional space in which to represent the graph.
    References
    ----------
    .. [1] Sussman, D.L., Tang, M., Fishkind, D.E., Priebe, C.E.  "A
       Consistent Adjacency Spectral Embedding for Stochastic Blockmodel Graphs,"
       Journal of the American Statistical Association, Vol. 107(499), 2012
    .. [2] Zhu, M. and Ghodsi, A. (2006).
        Automatic dimensionality selection from the scree plot via the use of
        profile likelihood. Computational Statistics & Data Analysis, 51(2),
        pp.918-930.
    .. [3] Levin, K., Roosta-Khorasani, F., Mahoney, M. W., & Priebe, C. E. (2018).
        Out-of-sample extension of graph adjacency spectral embedding. PMLR: Proceedings
        of Machine Learning Research, 80, 2975-2984.
    """

    def __init__(
        self,
        n_components=None,
        n_elbows=2,
        algorithm="ASE",
        svd_solver="randomized",
        check_lcc=True,
        regularizer=True,
        concat=False,
    ):

        self.n_components = n_components
        self.n_elbows = n_elbows
        self.algorithm = algorithm
        self.svd_solver = svd_solver
        self.check_lcc = check_lcc
        self.regularizer = regularizer
        self.concat = concat

    def fit(self, X, y=None):
        """
        Fit model to input graph

        Parameters
        ----------
        X : {array-like, sparse matrix} shape (n_samples, n_samples)
            Adjacency matrix of graph to embed.

        y: Ignored

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        A = self._validate_data(
            X, accept_sparse="csr", ensure_min_samples=2, estimator=self
        )
        if A.shape[0] != A.shape[1]:
            raise ValueError("Matrix input must be square")

        if self.check_lcc:
            directed = not _is_symmetric(A)
            n_components = connected_components(
                A, directed, connection="weak", return_labels=False
            )
            if n_components != 1:
                msg = (
                    "Input graph is not fully connected. Results may not"
                    + "be optimal. You can compute the largest connected component by"
                    + "using ``graspologic.utils.largest_connected_component``."
                )
                warnings.warn(msg, UserWarning)

        self.n_features_in_ = A.shape[0]

        if isinstance(self.algorithm, str):
            if self.algorithm.lower() not in {"ase", "lse"}:
                raise ValueError(
                    ("%s is not a valid embedding method. Expected " "'ASE' or 'LSE'")
                    % self.algorithm
                )
        else:
            raise TypeError('"algorithm" must be of type string')

        # reduces the dimension of adjacency matrix using the desired embedding method.
        if self.algorithm.lower() == "ase" and self.regularizer:
            A = _augment_diagonal(A)
        elif self.algorithm.lower() == "lse":
            A = _to_laplacian(A, regularizer=self.regularizer)

        U, D, V = _selectSVD(
            A,
            n_components=self.n_components,
            n_elbows=self.n_elbows,
            svd_solver=self.svd_solver,
        )

        self.n_components_ = D.size
        self.singular_values_ = D
        self.latent_left_ = U @ np.diag(np.sqrt(D))
        if not _is_symmetric(A):
            self.latent_right_ = V.T @ np.diag(np.sqrt(D))
        else:
            self.latent_right_ = None

        # for out-of-sample
        inv_eigs = np.diag(1 / self.singular_values_)
        self._pinv_left = self.latent_left_ @ inv_eigs
        if self.latent_right_ is not None:
            self._pinv_right = self.latent_right_ @ inv_eigs

        self.is_fitted_ = True
        return self

    def fit_transform(self, X, y=None):
        """
        Fit model to input graph

        Parameters
        ----------
        X : {array-like, sparse matrix} shape (n_samples, n_samples)
            Adjacency matrix of graph to embed.

        y: Ignored

        Returns
        -------
        X_new : array-like of shape (n_samples, n_components)
        """
        self.fit(X)
        if self.latent_right_ is None:
            return self.latent_left_
        else:
            if self.concat:
                return np.concatenate((self.latent_left_, self.latent_right_), axis=1)
            else:
                return self.latent_left_, self.latent_right_

    def transform(self, X):
        """
        Obtain latent positions from an adjacency matrix or matrix of out-of-sample
        vertices. For more details on transforming out-of-sample vertices, see the
        :ref:`tutorials <embed_tutorials>`. For mathematical background, see [2].

        Parameters
        ----------
        X : array-like or tuple, original shape or (n_oos_vertices, n_vertices).
            The original fitted matrix ("graph" in fit) or new out-of-sample data.
            If ``X`` is the original fitted matrix, returns a matrix close to
            ``self.fit_transform(X)``.

            If ``X`` is an out-of-sample matrix, n_oos_vertices is the number
            of new vertices, and n_vertices is the number of vertices in the
            original graph. If tuple, graph is directed and ``X[0]`` contains
            edges from out-of-sample vertices to in-sample vertices.

        Returns
        -------
        array_like or tuple, shape (n_oos_vertices, n_components)
        or (n_vertices, n_components).
            Array of latent positions. Transforms the fitted matrix if it was passed
            in.

            If ``X`` is an array or tuple containing adjacency vectors corresponding to
            new nodes, returns the estimated latent positions for the new out-of-sample
            adjacency vectors.

            If undirected, returns array.
            If directed, returns ``(X_out, X_in)``, where ``X_out`` contains
            latent positions corresponding to nodes with edges from out-of-sample
            vertices to in-sample vertices.

        Notes
        -----
        If the matrix was diagonally augmented (e.g., ``self.regularizer`` was True),
        ``fit`` followed by ``transform`` will produce a slightly different matrix than
        ``fit_transform``.

        To get the original embedding, using ``fit_transform`` is recommended. In the
        directed case, if A is the original in-sample adjacency matrix, the tuple
        (A.T, A) will need to be passed to ``transform`` if you do not wish to use
        ``fit_transform``.
        """

        # checks
        check_is_fitted(self, "is_fitted_")
        directed = self.latent_right_ is not None

        # correct types?
        if directed and not isinstance(X, tuple):
            if X.shape[0] == X.shape[1]:  # in case original matrix was passed
                msg = """A square matrix A was passed to ``transform`` in the directed case.
                If this was the original in-sample matrix, either use ``fit_transform``
                or pass a tuple (A.T, A). If this was an out-of-sample matrix, directed
                graphs require a tuple (X_out, X_in)."""
                raise TypeError(msg)
            else:
                msg = """Directed graphs require a tuple (X_out, X_in) for
                out-of-sample transforms."""
                raise TypeError(msg)
        if not directed and not isinstance(X, np.ndarray):
            raise TypeError("Undirected graphs require array input")

        # correct shape in y?
        latent_rows = self.latent_left_.shape[0]
        _X = X[0] if directed else X
        X_cols = _X.shape[-1]
        if _X.ndim > 2:
            raise ValueError("out-of-sample vertex must be 1d or 2d")
        if latent_rows != X_cols:
            msg = "out-of-sample vertex must be shape (n_oos_vertices, n_vertices)"
            raise ValueError(msg)

        # workhorse code
        if not directed:
            return X @ self._pinv_left
        elif directed:
            return X[1] @ self._pinv_right, X[0] @ self._pinv_left
