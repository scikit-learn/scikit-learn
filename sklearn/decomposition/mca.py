""" Multiple Correspondance Analysis
not possible without Vivek Yadak's help
http://vxy10.github.io/2016/06/10/intro-MCA/
"""

# TODO: Handle Sparse Arrays?
# TODO: handle benzecri correction
# TODO: Handle greenacre correction
# TODO: derive n_components if not provided?
# TODO: Add documentation
# TODO: Add lots of comments
# TODO: Add BaseClass
# TODO: add .score methods
# TODO: Make Expl_var accessible
# TODO: SVD svd_solver parameter

from scipy.linalg import diagsvd
from numpy.linalg import svd
import numpy as np

import functools
def Matrix_mult(*args):
    """An internal method to multiply matrices."""
    return functools.reduce(np.dot, args)

class MCA():
    """ Multiple Correspondance Analysis
    """

    def __init__(self, n_components=None):
        self.n_components = n_components

    def fit(self, X, y=None):
        """
        X : array-like, shape (n_samples, n_features)
        """
        self._fit(X)
        return self

    def _fit(self, X, y=None):
        # determine if we can call full _fit
        # if smaller than 500 dimensions

        return self._fit_full(X, self.n_components)


    def _fit_full(self, X, n_components):
        n_samples, n_features = X.shape
        i_sup = X.tail(1)

        N_all = np.sum(X)
        Z = X / N_all

        # Get 2 vectors corresponding to the sum of rows and colums.
        Sum_r = np.sum(Z,axis=1)
        Sum_c = np.sum(Z,axis=0)

        # Compute residual matrix by subtracting the expected indicator matrix
        # (outer product of rows and columns sums computed in step 2)
        Z_expected = np.outer(Sum_r, Sum_c)
        Z_residual = Z - Z_expected

        # Scale residual by the square root of column and row sums.
        D_r = np.diag(Sum_r)
        D_c = np.diag(Sum_c)
        D_r_sqrt_mi = np.sqrt(np.diag(Sum_r**-1))
        D_c_sqrt_mi = np.sqrt(np.diag(Sum_c**-1))

        # Apply SVD to the residual matrix to compute left- (row) and right- (column) eigenvectors.
        MCA_mat = Matrix_mult(D_r_sqrt_mi,Z_residual,D_c_sqrt_mi)
        P,S,Q = svd(MCA_mat) # most costly part, whole matrices set to false
        #S_d = diagsvd(S, 6, 22)
        S_d = diagsvd(S, n_samples, n_features)
        sum_mca = np.sum((Matrix_mult(P,S_d,Q)-MCA_mat)**2)
        # print 'Difference between SVD and the MCA matrix is %0.2f' % sum_mca

        F = Matrix_mult(D_r_sqrt_mi,P,S_d) ## Column Space, contains linear combinations of columns
        G = Matrix_mult(D_c_sqrt_mi,Q.T,S_d.T) ## Row space, contains linear combinations of rows
        Lam = S**2
        Expl_var = Lam/np.sum(Lam)

        # print 'Eigen values are ', Lam
        # print 'Explained variance of eigen vectors are ', Expl_var
        # Choose number of dimensions

        # Benzecri correction
        # TODO Check for this or the other correction
        K = 10
        E = np.array([(K/(K-1.)*(lm - 1./K))**2 if lm > 1./K else 0 for lm in S**2])
        Expl_var_bn = E/np.sum(E)

        dim = self.n_components
        d_c = np.linalg.norm(G, axis=1)**2 # Same as np.diag(G*G.T)
        d_r = np.linalg.norm(F, axis=1)**2 # Same as np.diag(F*F.T)
        # Cosine distance between factors and first 4 components
        CosDist_c = np.apply_along_axis(lambda x: x/d_c, 0, G[:, :dim]**2)
        # Cosine distance between factors and first dim components
        CosDist_r = np.apply_along_axis(lambda x: x/d_r, 0, F[:, :dim]**2)
        # Cosine distance between factors and first dim components
        Cont_c = np.apply_along_axis(lambda x: x, 0, G[:, :dim]**2)
        # Cosine distance between factors and first dim components
        Cont_r = np.apply_along_axis(lambda x: x, 0, F[:, :dim]**2)

        X_pjn = []
        for i in np.arange(0,6):
            X_pjn.append(np.dot(X.iloc[i],G[:,:dim])/S[:dim]/10)

        X_pjn = np.asarray(X_pjn)

        # THE RESULT
        X_pjn = Matrix_mult(X.values,G[:,:dim])/S[:dim]/10
        X_i = -Matrix_mult(i_sup,G[:,:dim])/S[:dim]/10
        return X_pjn

    def fit_transform(self, X, y=None):
        X = self._fit(X)
        return X
