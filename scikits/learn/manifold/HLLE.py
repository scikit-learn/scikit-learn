"""Locally Linear Embedding"""

# Author: Jake Vanderplas <vanderplas@astro.washington.edu>

import numpy as np
from ..base import BaseEstimator
from ..neighbors import kneighbors_graph
from ..ball_tree import BallTree


def hessian_locally_linear_embedding(
    X, n_neighbors, out_dim, reg=1e-3, eigen_solver='lobpcg', tol=1e-6,
    max_iter=100):
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

    Notes
    -----
      This follows the Hessian Eigenmapping algorithm outlined in [1].  
      Hessian Eigenmapping uses hessian weights in each neighborhood
      to recover a much more robust projection than LLE.  The additional
      constraint on the number of neighbors (k>Ndims), as well as the
      large computational cost, make HLLE unsuited for high-dimensional
      problems.

    References
    ----------
      [1] Donoho, D & Grimes, C. Hessian eigenmaps: Locally linear embedding 
          techniques for high-dimensional data. Proc Natl Acad Sci U S A. 
          100:5591 (2003)
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
    dp = d_out*(d_out+1)/2
        
    X = np.asarray(X)
    
    assert d_out < d_in
    assert k > d_out + dp #note that dp = d_out*(d_out+1)/2
    assert k < N

    W = np.zeros((N,N*dp),dtype=np.float)

    Yi = np.ones((k,1+d_out+dp),dtype=np.float)

    for i in range(N):
        Gi = X[neighbors[i]]
        Gi -= Gi.mean(0)

        U,sig,VT = np.linalg.svd(Gi,full_matrices=0)

        Yi[:,1:1+d_out] = U[:,:d_out]

        #build hessian estimator
        ct = 1+d_out
        for mm in range(d_out):
            for nn in range(mm,d_out):
                Yi[:,ct] = U[:,mm]*U[:,nn]
                ct+=1

        Q,R = np.linalg.qr(Yi)
           
        w = Q[:,d_out+1:]
        S = w.sum(0)

        S[np.where(abs(S)<1E-4)] = 1
        W[neighbors[i],i*dp:(i+1)*dp] = w/S

    if eigen_solver == 'dense':
        if True:
            eigen_vectors,sig,VT = np.linalg.svd(W,full_matrices=0)
            eigen_values = sig**2
            index = np.argsort(sig)[1:d_out+1]
        else:
            import scipy.linalg
            C = np.dot(W,W.T)
            eigen_values, eigen_vectors = scipy.linalg.eigh(
                C, eigvals=(1, out_dim + 1), overwrite_a=True)
            index = np.argsort(np.abs(eigen_values))
        Y = eigen_vectors[:,index]
    
    elif eigen_solver == 'lobpcg':
        from scipy.sparse import linalg, eye, csr_matrix
        C = np.dot(W,W.T)
        
        C = csr_matrix(C)
    
        # initial approximation to the eigenvectors
        X = np.random.rand(C.shape[0], out_dim+1)
        try:
            ml = pyamg.smoothed_aggregation_solver(C, symmetry='symmetric')
        except TypeError:
            ml = pyamg.smoothed_aggregation_solver(C, mat_flag='symmetric')
        prec = ml.aspreconditioner()

        # compute eigenvalues and eigenvectors with LOBPCG
        eigen_values, eigen_vectors = linalg.lobpcg(
            C, X, M=prec, largest=False, tol=tol, maxiter=max_iter)

        index = np.argsort(eigen_values)
        Y = eigen_vectors[:, index[1:]]
        eigen_values = eigen_values[index[1:]]

    else:
        raise NotImplementedError('Method %s not implemented' % eigen_solver)

    Y *= np.sqrt(N)

    #-----------------------------------------------
    # Normalize Y
    #  we need R = (Y.T*Y)^(-1/2)
    #   do this with an SVD of Y
    #      Y = U*sig*V.T
    #      Y.T*Y = (V*sig.T*U.T) * (U*sig*V.T)
    #            = U*(sig*sig.T)*U.T
    #   so
    #      R = V * sig^-1 * V.T
    #-----------------------------------------------
    #U,sig,VT = numpy.linalg.svd(Y,full_matrices=0)
    #del U
    #S = numpy.matrix(numpy.diag(sig**-1))
    #R = VT.T * S * VT
    #return numpy.array(Y*R)

    U,sig,VT = np.linalg.svd(Y.T,full_matrices=0)
    S = np.diag(sig**-1)
    R = np.dot(VT.T, np.dot(S,VT))
    return np.dot(Y.T,R).T, np.sum(eigen_vectors)



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
