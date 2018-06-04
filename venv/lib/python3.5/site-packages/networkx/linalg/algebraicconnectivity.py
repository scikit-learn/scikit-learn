# -*- coding: utf-8 -*-
# Copyright (C) 2014 ysitu <ysitu@users.noreply.github.com>
# All rights reserved.
# BSD license.
#
# Author: ysitu <ysitu@users.noreply.github.com>
"""
Algebraic connectivity and Fiedler vectors of undirected graphs.
"""
from functools import partial
import networkx as nx
from networkx.utils import not_implemented_for
from networkx.utils import reverse_cuthill_mckee_ordering

try:
    from numpy import array, asmatrix, asarray, dot, ndarray, ones, sqrt, zeros
    from numpy.linalg import norm, qr
    from numpy.random import normal
    from scipy.linalg import eigh, inv
    from scipy.sparse import csc_matrix, spdiags
    from scipy.sparse.linalg import eigsh, lobpcg
    __all__ = ['algebraic_connectivity', 'fiedler_vector', 'spectral_ordering']
except ImportError:
    __all__ = []

try:
    from scipy.linalg.blas import dasum, daxpy, ddot
except ImportError:
    if __all__:
        # Make sure the imports succeeded.
        # Use minimal replacements if BLAS is unavailable from SciPy.
        dasum = partial(norm, ord=1)
        ddot = dot

        def daxpy(x, y, a):
            y += a * x
            return y


class _PCGSolver(object):
    """Preconditioned conjugate gradient method.

    To solve Ax = b:
        M = A.diagonal() # or some other preconditioner
        solver = _PCGSolver(lambda x: A * x, lambda x: M * x)
        x = solver.solve(b)

    The inputs A and M are functions which compute
    matrix multiplication on the argument.
    A - multiply by the matrix A in Ax=b
    M - multiply by M, the preconditioner surragate for A

    Warning: There is no limit on number of iterations.
    """

    def __init__(self, A, M):
        self._A = A
        self._M = M or (lambda x: x.copy())

    def solve(self, B, tol):
        B = asarray(B)
        X = ndarray(B.shape, order='F')
        for j in range(B.shape[1]):
            X[:, j] = self._solve(B[:, j], tol)
        return X

    def _solve(self, b, tol):
        A = self._A
        M = self._M
        tol *= dasum(b)
        # Initialize.
        x = zeros(b.shape)
        r = b.copy()
        z = M(r)
        rz = ddot(r, z)
        p = z.copy()
        # Iterate.
        while True:
            Ap = A(p)
            alpha = rz / ddot(p, Ap)
            x = daxpy(p, x, a=alpha)
            r = daxpy(Ap, r, a=-alpha)
            if dasum(r) < tol:
                return x
            z = M(r)
            beta = ddot(r, z)
            beta, rz = beta / rz, beta
            p = daxpy(p, z, a=beta)


class _CholeskySolver(object):
    """Cholesky factorization.

    To solve Ax = b:
        solver = _CholeskySolver(A)
        x = solver.solve(b)

    optional argument `tol` on solve method is ignored but included
    to match _PCGsolver API.
    """

    def __init__(self, A):
        if not self._cholesky:
            raise nx.NetworkXError('Cholesky solver unavailable.')
        self._chol = self._cholesky(A)

    def solve(self, B, tol=None):
        return self._chol(B)

    try:
        from scikits.sparse.cholmod import cholesky
        _cholesky = cholesky
    except ImportError:
        _cholesky = None


class _LUSolver(object):
    """LU factorization.

    To solve Ax = b:
        solver = _LUSolver(A)
        x = solver.solve(b)

    optional argument `tol` on solve method is ignored but included
    to match _PCGsolver API.
    """

    def __init__(self, A):
        if not self._splu:
            raise nx.NetworkXError('LU solver unavailable.')
        self._LU = self._splu(A)

    def solve(self, B, tol=None):
        B = asarray(B)
        X = ndarray(B.shape, order='F')
        for j in range(B.shape[1]):
            X[:, j] = self._LU.solve(B[:, j])
        return X

    try:
        from scipy.sparse.linalg import splu
        _splu = partial(splu, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.,
                        options={'Equil': True, 'SymmetricMode': True})
    except ImportError:
        _splu = None


def _preprocess_graph(G, weight):
    """Compute edge weights and eliminate zero-weight edges.
    """
    if G.is_directed():
        H = nx.MultiGraph()
        H.add_nodes_from(G)
        H.add_weighted_edges_from(((u, v, e.get(weight, 1.))
                                   for u, v, e in G.edges(data=True)
                                   if u != v), weight=weight)
        G = H
    if not G.is_multigraph():
        edges = ((u, v, abs(e.get(weight, 1.)))
                 for u, v, e in G.edges(data=True) if u != v)
    else:
        edges = ((u, v, sum(abs(e.get(weight, 1.)) for e in G[u][v].values()))
                 for u, v in G.edges() if u != v)
    H = nx.Graph()
    H.add_nodes_from(G)
    H.add_weighted_edges_from((u, v, e) for u, v, e in edges if e != 0)
    return H


def _rcm_estimate(G, nodelist):
    """Estimate the Fiedler vector using the reverse Cuthill-McKee ordering.
    """
    G = G.subgraph(nodelist)
    order = reverse_cuthill_mckee_ordering(G)
    n = len(nodelist)
    index = dict(zip(nodelist, range(n)))
    x = ndarray(n, dtype=float)
    for i, u in enumerate(order):
        x[index[u]] = i
    x -= (n - 1) / 2.
    return x


def _tracemin_fiedler(L, X, normalized, tol, method):
    """Compute the Fiedler vector of L using the TraceMIN-Fiedler algorithm.

    The Fiedler vector of a connected undirected graph is the eigenvector
    corresponding to the second smallest eigenvalue of the Laplacian matrix of
    of the graph. This function starts with the Laplacian L, not the Graph.

    Parameters
    ----------
    L : Laplacian of a possibly weighted or normalized, but undirected graph

    X : Initial guess for a solution. Usually a matrix of random numbers.
        This function allows more than one column in X to identify more than
        one eigenvector if desired.

    normalized : bool
        Whether the normalized Laplacian matrix is used.

    tol : float
        Tolerance of relative residual in eigenvalue computation.
        Warning: There is no limit on number of iterations.

    method : string
        Should be 'tracemin_pcg', 'tracemin_chol' or 'tracemin_lu'.
        Otherwise exception is raised.

    Returns
    -------
    sigma, X : Two NumPy arrays of floats.
        The lowest eigenvalues and corresponding eigenvectors of L.
        The size of input X determines the size of these outputs.
        As this is for Fiedler vectors, the zero eigenvalue (and
        constant eigenvector) are avoided.
    """
    n = X.shape[0]

    if normalized:
        # Form the normalized Laplacian matrix and determine the eigenvector of
        # its nullspace.
        e = sqrt(L.diagonal())
        D = spdiags(1. / e, [0], n, n, format='csr')
        L = D * L * D
        e *= 1. / norm(e, 2)

    if normalized:
        def project(X):
            """Make X orthogonal to the nullspace of L.
            """
            X = asarray(X)
            for j in range(X.shape[1]):
                X[:, j] -= dot(X[:, j], e) * e
    else:
        def project(X):
            """Make X orthogonal to the nullspace of L.
            """
            X = asarray(X)
            for j in range(X.shape[1]):
                X[:, j] -= X[:, j].sum() / n

    if method == 'tracemin_pcg':
        D = L.diagonal().astype(float)
        solver = _PCGSolver(lambda x: L * x, lambda x: D * x)
    elif method == 'tracemin_chol' or method == 'tracemin_lu':
        # Convert A to CSC to suppress SparseEfficiencyWarning.
        A = csc_matrix(L, dtype=float, copy=True)
        # Force A to be nonsingular. Since A is the Laplacian matrix of a
        # connected graph, its rank deficiency is one, and thus one diagonal
        # element needs to modified. Changing to infinity forces a zero in the
        # corresponding element in the solution.
        i = (A.indptr[1:] - A.indptr[:-1]).argmax()
        A[i, i] = float('inf')
        if method == 'tracemin_chol':
            solver = _CholeskySolver(A)
        else:
            solver = _LUSolver(A)
    else:
        raise nx.NetworkXError('Unknown linear system solver: ' + method)

    # Initialize.
    Lnorm = abs(L).sum(axis=1).flatten().max()
    project(X)
    W = asmatrix(ndarray(X.shape, order='F'))

    while True:
        # Orthonormalize X.
        X = qr(X)[0]
        # Compute interation matrix H.
        W[:, :] = L * X
        H = X.T * W
        sigma, Y = eigh(H, overwrite_a=True)
        # Compute the Ritz vectors.
        X *= Y
        # Test for convergence exploiting the fact that L * X == W * Y.
        res = dasum(W * asmatrix(Y)[:, 0] - sigma[0] * X[:, 0]) / Lnorm
        if res < tol:
            break
        # Compute X = L \ X / (X' * (L \ X)).
        # L \ X can have an arbitrary projection on the nullspace of L,
        # which will be eliminated.
        W[:, :] = solver.solve(X, tol)
        X = (inv(W.T * X) * W.T).T  # Preserves Fortran storage order.
        project(X)

    return sigma, asarray(X)


def _get_fiedler_func(method):
    """Return a function that solves the Fiedler eigenvalue problem.
    """
    if method == "tracemin":  # old style keyword <v2.1
        method = "tracemin_pcg"
    if method in ("tracemin_pcg", "tracemin_chol", "tracemin_lu"):
        def find_fiedler(L, x, normalized, tol):
            q = 1 if method == 'tracemin_pcg' else min(4, L.shape[0] - 1)
            X = asmatrix(normal(size=(q, L.shape[0]))).T
            sigma, X = _tracemin_fiedler(L, X, normalized, tol, method)
            return sigma[0], X[:, 0]
    elif method == 'lanczos' or method == 'lobpcg':
        def find_fiedler(L, x, normalized, tol):
            L = csc_matrix(L, dtype=float)
            n = L.shape[0]
            if normalized:
                D = spdiags(1. / sqrt(L.diagonal()), [0], n, n, format='csc')
                L = D * L * D
            if method == 'lanczos' or n < 10:
                # Avoid LOBPCG when n < 10 due to
                # https://github.com/scipy/scipy/issues/3592
                # https://github.com/scipy/scipy/pull/3594
                sigma, X = eigsh(L, 2, which='SM', tol=tol,
                                 return_eigenvectors=True)
                return sigma[1], X[:, 1]
            else:
                X = asarray(asmatrix(x).T)
                M = spdiags(1. / L.diagonal(), [0], n, n)
                Y = ones(n)
                if normalized:
                    Y /= D.diagonal()
                sigma, X = lobpcg(L, X, M=M, Y=asmatrix(Y).T, tol=tol,
                                  maxiter=n, largest=False)
                return sigma[0], X[:, 0]
    else:
        raise nx.NetworkXError("unknown method '%s'." % method)

    return find_fiedler


@not_implemented_for('directed')
def algebraic_connectivity(G, weight='weight', normalized=False, tol=1e-8,
                           method='tracemin_pcg'):
    """Return the algebraic connectivity of an undirected graph.

    The algebraic connectivity of a connected undirected graph is the second
    smallest eigenvalue of its Laplacian matrix.

    Parameters
    ----------
    G : NetworkX graph
        An undirected graph.

    weight : object, optional (default: None)
        The data key used to determine the weight of each edge. If None, then
        each edge has unit weight.

    normalized : bool, optional (default: False)
        Whether the normalized Laplacian matrix is used.

    tol : float, optional (default: 1e-8)
        Tolerance of relative residual in eigenvalue computation.

    method : string, optional (default: 'tracemin_pcg')
        Method of eigenvalue computation. It must be one of the tracemin
        options shown below (TraceMIN), 'lanczos' (Lanczos iteration)
        or 'lobpcg' (LOBPCG).

        The TraceMIN algorithm uses a linear system solver. The following
        values allow specifying the solver to be used.

        =============== ========================================
        Value           Solver
        =============== ========================================
        'tracemin_pcg'  Preconditioned conjugate gradient method
        'tracemin_chol' Cholesky factorization
        'tracemin_lu'   LU factorization
        =============== ========================================

    Returns
    -------
    algebraic_connectivity : float
        Algebraic connectivity.

    Raises
    ------
    NetworkXNotImplemented
        If G is directed.

    NetworkXError
        If G has less than two nodes.

    Notes
    -----
    Edge weights are interpreted by their absolute values. For MultiGraph's,
    weights of parallel edges are summed. Zero-weighted edges are ignored.

    To use Cholesky factorization in the TraceMIN algorithm, the
    :samp:`scikits.sparse` package must be installed.

    See Also
    --------
    laplacian_matrix
    """
    if len(G) < 2:
        raise nx.NetworkXError('graph has less than two nodes.')
    G = _preprocess_graph(G, weight)
    if not nx.is_connected(G):
        return 0.

    L = nx.laplacian_matrix(G)
    if L.shape[0] == 2:
        return 2. * L[0, 0] if not normalized else 2.

    find_fiedler = _get_fiedler_func(method)
    x = None if method != 'lobpcg' else _rcm_estimate(G, G)
    sigma, fiedler = find_fiedler(L, x, normalized, tol)
    return sigma


@not_implemented_for('directed')
def fiedler_vector(G, weight='weight', normalized=False, tol=1e-8,
                   method='tracemin_pcg'):
    """Return the Fiedler vector of a connected undirected graph.

    The Fiedler vector of a connected undirected graph is the eigenvector
    corresponding to the second smallest eigenvalue of the Laplacian matrix of
    of the graph.

    Parameters
    ----------
    G : NetworkX graph
        An undirected graph.

    weight : object, optional (default: None)
        The data key used to determine the weight of each edge. If None, then
        each edge has unit weight.

    normalized : bool, optional (default: False)
        Whether the normalized Laplacian matrix is used.

    tol : float, optional (default: 1e-8)
        Tolerance of relative residual in eigenvalue computation.

    method : string, optional (default: 'tracemin_pcg')
        Method of eigenvalue computation. It must be one of the tracemin
        options shown below (TraceMIN), 'lanczos' (Lanczos iteration)
        or 'lobpcg' (LOBPCG).

        The TraceMIN algorithm uses a linear system solver. The following
        values allow specifying the solver to be used.

        =============== ========================================
        Value           Solver
        =============== ========================================
        'tracemin_pcg'  Preconditioned conjugate gradient method
        'tracemin_chol' Cholesky factorization
        'tracemin_lu'   LU factorization
        =============== ========================================

    Returns
    -------
    fiedler_vector : NumPy array of floats.
        Fiedler vector.

    Raises
    ------
    NetworkXNotImplemented
        If G is directed.

    NetworkXError
        If G has less than two nodes or is not connected.

    Notes
    -----
    Edge weights are interpreted by their absolute values. For MultiGraph's,
    weights of parallel edges are summed. Zero-weighted edges are ignored.

    To use Cholesky factorization in the TraceMIN algorithm, the
    :samp:`scikits.sparse` package must be installed.

    See Also
    --------
    laplacian_matrix
    """
    if len(G) < 2:
        raise nx.NetworkXError('graph has less than two nodes.')
    G = _preprocess_graph(G, weight)
    if not nx.is_connected(G):
        raise nx.NetworkXError('graph is not connected.')

    if len(G) == 2:
        return array([1., -1.])

    find_fiedler = _get_fiedler_func(method)
    L = nx.laplacian_matrix(G)
    x = None if method != 'lobpcg' else _rcm_estimate(G, G)
    sigma, fiedler = find_fiedler(L, x, normalized, tol)
    return fiedler


def spectral_ordering(G, weight='weight', normalized=False, tol=1e-8,
                      method='tracemin_pcg'):
    """Compute the spectral_ordering of a graph.

    The spectral ordering of a graph is an ordering of its nodes where nodes
    in the same weakly connected components appear contiguous and ordered by
    their corresponding elements in the Fiedler vector of the component.

    Parameters
    ----------
    G : NetworkX graph
        A graph.

    weight : object, optional (default: None)
        The data key used to determine the weight of each edge. If None, then
        each edge has unit weight.

    normalized : bool, optional (default: False)
        Whether the normalized Laplacian matrix is used.

    tol : float, optional (default: 1e-8)
        Tolerance of relative residual in eigenvalue computation.

    method : string, optional (default: 'tracemin_pcg')
        Method of eigenvalue computation. It must be one of the tracemin
        options shown below (TraceMIN), 'lanczos' (Lanczos iteration)
        or 'lobpcg' (LOBPCG).

        The TraceMIN algorithm uses a linear system solver. The following
        values allow specifying the solver to be used.

        =============== ========================================
        Value           Solver
        =============== ========================================
        'tracemin_pcg'  Preconditioned conjugate gradient method
        'tracemin_chol' Cholesky factorization
        'tracemin_lu'   LU factorization
        =============== ========================================

    Returns
    -------
    spectral_ordering : NumPy array of floats.
        Spectral ordering of nodes.

    Raises
    ------
    NetworkXError
        If G is empty.

    Notes
    -----
    Edge weights are interpreted by their absolute values. For MultiGraph's,
    weights of parallel edges are summed. Zero-weighted edges are ignored.

    To use Cholesky factorization in the TraceMIN algorithm, the
    :samp:`scikits.sparse` package must be installed.

    See Also
    --------
    laplacian_matrix
    """
    if len(G) == 0:
        raise nx.NetworkXError('graph is empty.')
    G = _preprocess_graph(G, weight)

    find_fiedler = _get_fiedler_func(method)
    order = []
    for component in nx.connected_components(G):
        size = len(component)
        if size > 2:
            L = nx.laplacian_matrix(G, component)
            x = None if method != 'lobpcg' else _rcm_estimate(G, component)
            sigma, fiedler = find_fiedler(L, x, normalized, tol)
            sort_info = zip(fiedler, range(size), component)
            order.extend(u for x, c, u in sorted(sort_info))
        else:
            order.extend(component)

    return order


# fixture for nose tests
def setup_module(module):
    from nose import SkipTest
    try:
        import numpy
        import scipy.sparse
    except ImportError:
        raise SkipTest('SciPy not available.')
