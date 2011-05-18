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
        
    X = np.asarray(X).T

    assert d_out < d_in
    assert k >= d_out
    assert k < N

    #some variables to hold needed values
    rho = np.zeros(N)
    w_reg = np.zeros([N,k])
    evals = np.zeros([N,k])
    V = np.zeros((N,k,k))

    Phi = np.zeros( (N,N) )

    #some functions to simplify the code
    one = lambda d: np.ones((d,1))

    for i in range(N):
        #find regularized weights: this is like normal LLE
        Gi = X[:,neighbors[i]] - X[:,i:i+1]
        Qi = np.dot(Gi.T,Gi)

        Qi.flat[::k+1] += 1E-3 * Qi.trace()

        w_reg[i] = np.linalg.solve(Qi,np.ones(k))
        w_reg[i] /= w_reg[i].sum()

        #find the eigenvectors and eigenvalues of Gi.T*Gi
        # using SVD
        # we want V[i] to be a [k x k] matrix, where the columns
        # are eigenvectors of Gi^T * G
        V[i],sig,UT = np.linalg.svd(Gi.T)
        evals[i][:len(sig)] = sig**2

        #compute rho_i : this is used to calculate eta, the
        # cutoff used to determine the size of the "almost null"
        # space of the local covariance matrices.
        rho[i] = (evals[i,d_out:]).sum() / (evals[i,:d_out]).sum()

    #find eta - the median of the N rho values
    rho.sort()
    eta = rho[int(N/2)]

    #The next loop calculates Phi.
    # This is the [N x N] matrix whose null space is the desired embedding
    for i in range(N):
        #determine si - the size of the largest set of eigenvalues
        # of Qi such that satisfies:
        #    sum(in_set)/sum(not_in_set) < eta
        # with the constraint that 0<si<=k-d_out

        si = 1
        while si < k-d_out:
            this_eta = sum( evals[i,k-si:] ) / sum( evals[i,:k-si] )
            if this_eta > eta:
                if(si!=1): si -= 1
                break
            else:
                si+=1

        #select bottom si eigenvectors of Qi
        # and calculate alpha
        Vi = V[i][:,k-si:]
        alpha_i = np.linalg.norm( Vi.sum(0) )/np.sqrt(si)

        #compute Householder matrix which satisfies
        #  Hi*Vi.T*ones(k) = alpha_i*ones(s)
        # using prescription from paper
        h = alpha_i * np.ones(si)[:,None] - np.dot(Vi.T,np.ones(k)[:,None])

        nh = np.linalg.norm(h)
        if nh < MLLE_tol:
            h = np.zeros( (si,1) )
        else:
            h /= nh

        Hi = np.identity(si) - 2*np.dot(h,h.T)

        Wi = np.dot(Vi,Hi) + (1-alpha_i) * w_reg[i][:,None]

        W_hat = np.zeros( (N,si) )
        W_hat[neighbors[i],:] = Wi
        W_hat[i]-=1
            
        Phi += np.dot(W_hat,W_hat.T)
    
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
