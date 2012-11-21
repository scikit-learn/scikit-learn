"""Temporary file to collate all useful functions
"""

# Covariance estimation methods
# -----------------------------
# input: samples x features
# returns samples x samples


# http://projecteuclid.org/DPubS?verb=Display&version=1.0&service=UI&handle=euclid.aoas/1267453933&page=record

# http://www.springerlink.com/content/w2g1j83430315572/

# http://strimmerlab.org/software/corpcor/


# Graph Laplacian creation
# ------------------------

# input: samples x samples
# returns samples x samples

# follow details from luxberg - tutorial on spectral clustering

# being implemented in ..utils.graph

# Eigenvalue estimation
# ---------------------

#input: samples x samples, n_eigs, smallest|biggest
#return: vectors x values


# from locally_linear
def null_space(M, k, k_skip=1, eigen_solver='arpack', tol=1E-6, max_iter=100,
               random_state=None):
    """
    Find the null space of a matrix M.

    Parameters
    ----------
    M : {array, matrix, sparse matrix, LinearOperator}
        Input covariance matrix: should be symmetric positive semi-definite

    k : integer
        Number of eigenvalues/vectors to return

    k_skip : integer, optional
        Number of low eigenvalues to skip.

    eigen_solver : string, {'auto', 'arpack', 'dense'}
        auto : algorithm will attempt to choose the best method for input data
        arpack : use arnoldi iteration in shift-invert mode.
                    For this method, M may be a dense matrix, sparse matrix,
                    or general linear operator.
        dense  : use standard dense matrix operations for the eigenvalue
                    decomposition.  For this method, M must be an array
                    or matrix type.  This method should be avoided for
                    large problems.

    tol : float, optional
        Tolerance for 'arpack' method.
        Not used if eigen_solver=='dense'.

    max_iter : maximum number of iterations for 'arpack' method
        not used if eigen_solver=='dense'

    random_state: numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.
        Used to determine the starting vector for arpack iterations
    """
    if eigen_solver == 'auto':
        if M.shape[0] > 200 and k + k_skip < 10:
            eigen_solver = 'arpack'
        else:
            eigen_solver = 'dense'

    if eigen_solver == 'arpack':
        random_state = check_random_state(random_state)
        v0 = random_state.rand(M.shape[0])
        try:
            eigen_values, eigen_vectors = eigsh(M, k + k_skip, sigma=0.0,
                                                tol=tol, maxiter=max_iter,
                                                v0=v0)
        except RuntimeError as msg:
            raise ValueError("Error in determining null-space with ARPACK. "
                             "Error message: '%s'. "
                             "Note that method='arpack' can fail when the "
                             "weight matrix is singular or otherwise "
                             "ill-behaved.  method='dense' is recommended. "
                             "See online documentation for more information."
            % msg)
        except:
            #let other errors pass through
            raise

        return eigen_vectors[:, k_skip:], np.sum(eigen_values[k_skip:])
    elif eigen_solver == 'dense':
        if hasattr(M, 'toarray'):
            M = M.toarray()
        eigen_values, eigen_vectors = eigh(
            M, eigvals=(k_skip, k + k_skip - 1), overwrite_a=True)
        index = np.argsort(np.abs(eigen_values))
        return eigen_vectors[:, index], np.sum(eigen_values)
    else:
        raise ValueError("Unrecognized eigen_solver '%s'" % eigen_solver)

# from ..cluster.spectral
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
        lambdas, diffusion_map = eigsh(-laplacian, k=n_components,
                                       sigma=1.0, which='LM')
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
