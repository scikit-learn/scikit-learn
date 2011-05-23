"""Locally Linear Embedding"""

# Author: Jake Vanderplas <vanderplas@astro.washington.edu>

import numpy as np
from ..base import BaseEstimator
from ..neighbors import kneighbors_graph
from ..ball_tree import BallTree


def modified_locally_linear_embedding(
    X, n_neighbors, out_dim, reg=1e-3, eigen_solver='lobpcg', tol=1e-6,
    max_iter=100, MLLE_tol=1E-12):
    """
    Perform a Modified Locally Linear Embedding (MLLE) analysis on the data.
    See reference [1]

    Parameters
    ----------
    X : array-like, shape [n_samples, n_features]

    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold.

    reg : float
        regularization constant, multiplies the trace of the local covariance
        matrix of the distances.

    eigen_solver : {'lobpcg', 'dense'}
        use the lobpcg eigensolver or a dense eigensolver based on LAPACK
        routines. The lobpcg solver is usually faster but depends on PyAMG.

    max_iter : integer
        maximum number of iterations for the lobpcg solver.

    Returns
    -------
    Y : array-like, shape [n_samples, out_dim]
        Embedding vectors.

    squared_error : float
       Reconstruction error for the embedding vectors. Equivalent to
       norm(Y - W Y, 'fro')**2, where W are the reconstruction weights.

    References
    ----------
      [1] Zhang,z & Wang, J. MLLE: Modified Locally Linear Embedding 
          Using Multiple Weights. 
          http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.70.382
    """

    if eigen_solver == 'lobpcg':
        try:
            import pyamg
        except ImportError:
            import warnings
            warnings.warn('amg was not found. Using (slow) dense eigensolver')
            eigen_solver = 'dense'

    balltree = BallTree(X)
    neighbors = balltree.query(X,k=n_neighbors+1,return_distance=False)
    neighbors = neighbors[:,1:]

    N,d_in = X.shape
    k = n_neighbors
    d_out = out_dim
        
    X = np.asarray(X)
    
    assert d_out < d_in
    assert k >= d_out
    assert k < N
    
    #find the eigenvectors and eigenvalues of each local covariance matrix
    # we want V[i] to be a [k x k] matrix, where the columns are eigenvectors
    V = np.zeros((N,k,k))
    nev = min(d_in,k)
    evals = np.zeros([N,nev])
    for i in range(N):
        V[i],evals[i],_ = np.linalg.svd(X[neighbors[i]] - X[i])
    evals**=2

    #find regularized weights: this is like normal LLE.
    # because we've already computed the SVD of each covariance matrix, 
    # it's faster to use this rather than np.linalg.solve
    reg = 1E-3 * evals.sum(1)

    tmp = np.dot(V.transpose(0,2,1),np.ones(k))
    tmp[:,:nev] /= evals+reg[:,None]
    tmp[:,nev:] /= reg[:,None]

    w_reg = np.zeros( (N,k) )
    for i in range(N):
        w_reg[i] = np.dot(V[i],tmp[i])
    w_reg /= w_reg.sum(1)[:,None]
    
    #calculate eta: the median of the ratio of small to large eigenvalues
    # across the points.  This is used to determine s_i, below
    rho = evals[:,d_out:].sum(1) / evals[:,:d_out].sum(1)
    eta = np.median(rho)

    #find s_i, the size of the "almost null space" for each point: 
    # this is the size of the largest set of eigenvalues
    # such that Sum[v; v in set]/Sum[v; v not in set] < eta
    s_range = np.zeros(N,dtype=int)
    evals_cumsum = np.cumsum(evals,1)
    eta_range    = evals_cumsum[:,-2:]/evals_cumsum[:,:-1] - 1
    for i in range(N):
        s_range[i] = k-2-np.searchsorted(eta_range[i],eta)

    #Now calculate Phi.
    # This is the [N x N] matrix whose null space is the desired embedding
    Phi = np.zeros( (N,N),dtype=np.float )
    for i in range(N):
        s_i = s_range[i]

        #select bottom s_i eigenvectors and calculate alpha
        Vi = V[i, :, k-s_i:]
        alpha_i = np.linalg.norm( Vi.sum(0) )/np.sqrt(s_i)

        #compute Householder matrix which satisfies
        #  Hi*Vi.T*ones(k) = alpha_i*ones(s)
        # using prescription from paper
        h = alpha_i * np.ones(s_i) - np.dot(Vi.T,np.ones(k))

        norm_h = np.linalg.norm(h)
        if norm_h < MLLE_tol: h *= 0
        else:                 h /= norm_h

        #Householder matrix is
        #  >> Hi = np.identity(s_i) - 2*np.outer(h,h)
        #Then the weight matrix is
        #  >> Wi = np.dot(Vi,Hi) + (1-alpha_i) * w_reg[i,:,None]
        #We do this much more efficiently:
        Wi = Vi - 2*np.outer(np.dot(Vi,h),h) + (1-alpha_i)*w_reg[i,:,None]

        #Update Phi as follows:
        # >> W_hat = np.zeros( (N,s_i) )
        # >> W_hat[neighbors[i],:] = Wi
        # >> W_hat[i] -= 1
        # >> Phi += np.dot(W_hat,W_hat.T)
        #We can do this much more efficiently:
        nbrs_x,nbrs_y = np.meshgrid(neighbors[i],neighbors[i])
        Phi[nbrs_x,nbrs_y]  += np.dot(Wi,Wi.T)
        Wi_sum1 = Wi.sum(1)
        Phi[i,neighbors[i]] -= Wi_sum1
        Phi[neighbors[i],i] -= Wi_sum1
        Phi[i,i]            += s_i
    
    if eigen_solver == 'dense':
        import scipy.linalg

        eigen_values, eigen_vectors = scipy.linalg.eigh(
            Phi, eigvals=(1, out_dim + 1), overwrite_a=True)
        index = np.argsort(np.abs(eigen_values))
        return eigen_vectors[:, index], np.sum(eigen_values)

    elif eigen_solver == 'lobpcg':
        from scipy.sparse import linalg, eye, csr_matrix
        
        Phi = csr_matrix(Phi)

        # initial approximation to the eigenvectors
        X = np.random.rand(Phi.shape[0], out_dim+1)
        try:
            ml = pyamg.smoothed_aggregation_solver(Phi, symmetry='symmetric')
        except TypeError:
            ml = pyamg.smoothed_aggregation_solver(Phi, mat_flag='symmetric')
        prec = ml.aspreconditioner()

        # compute eigenvalues and eigenvectors with LOBPCG
        eigen_values, eigen_vectors = linalg.lobpcg(
            Phi, X, M=prec, largest=False, tol=tol, maxiter=max_iter)

        index = np.argsort(eigen_values)
        return eigen_vectors[:, index[1:]], np.sum(eigen_values[index[1:]])

    else:
        raise NotImplementedError('Method %s not implemented' % eigen_solver)


class LocallyLinearEmbedding(BaseEstimator):
    """
    Locally Linear Embedding

    Parameters
    ----------
    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold

    Attributes
    ----------
    `embedding_vectors_` : array-like, shape [out_dim, n_samples]
        Stores the embedding vectors

    `reconstruction_error_` : float
        Reconstruction error associated with `embedding_vectors_`
    """

    def __init__(self, n_neighbors=5, out_dim=2):
        self.n_neighbors = n_neighbors
        self.out_dim = out_dim

    def fit(self, X, reg=1e-3, eigen_solver='lobpcg', tol=1e-6,
            max_iter=100, **params):
        """
        Compute the embedding vectors for data X.

        Parameters
        ----------
        X : array-like of shape [n_samples, n_features]
            training set.

        out_dim : integer
            number of coordinates for the manifold

        reg : float
            regularization constant, multiplies the trace of the local
            covariance matrix of the distances

        eigen_solver : {'lobpcg', 'dense'}
            use the lobpcg eigensolver or a dense eigensolver based on LAPACK
            routines. The lobpcg solver is usually faster but depends on PyAMG

        max_iter : integer
            maximum number of iterations for the lobpcg solver.

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)
        self.embedding_, self.reconstruction_error_ = \
            locally_linear_embedding(
                X, self.n_neighbors, self.out_dim, reg=reg,\
                eigen_solver=eigen_solver, tol=tol, max_iter=max_iter)
