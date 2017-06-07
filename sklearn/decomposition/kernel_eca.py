"""Kernel Entropy Components Analysis"""

# Author: Tobias Sterbak <sterbak-it@outlook.com>
# License: BSD 3 clause

import numpy as np
from scipy import linalg

from ..utils import check_random_state
from ..utils.arpack import eigsh
from ..utils.validation import check_is_fitted
from ..base import BaseEstimator, TransformerMixin
from ..preprocessing import KernelCenterer
from ..metrics.pairwise import pairwise_kernels

class KernelECA(BaseEstimator, TransformerMixin):
    """Kernel Entropy component analysis (KECA)

    Non-linear dimensionality reduction through the use of kernels (see
    :ref:`metrics`).

    Parameters
    ----------
    n_components: int or None
        Number of components. If None, all non-zero components are kept.

    kernel: "linear" | "poly" | "rbf" | "sigmoid" | "cosine" | "precomputed"
        Kernel.
        Default: "linear"

    degree : int, default=3
        Degree for poly kernels. Ignored by other kernels.

    gamma : float, optional
        Kernel coefficient for rbf and poly kernels. Default: 1/n_features.
        Ignored by other kernels.

    coef0 : float, optional
        Independent term in poly and sigmoid kernels.
        Ignored by other kernels.

    kernel_params : mapping of string to any, optional
        Parameters (keyword arguments) and values for kernel passed as
        callable object. Ignored by other kernels.


    eigen_solver: string ['auto'|'dense'|'arpack']
        Select eigensolver to use.  If n_components is much less than
        the number of training samples, arpack may be more efficient
        than the dense eigensolver.

    tol: float
        convergence tolerance for arpack.
        Default: 0 (optimal value will be chosen by arpack)

    max_iter : int
        maximum number of iterations for arpack
        Default: None (optimal value will be chosen by arpack)
		    
	random_state : int seed, RandomState instance, or None, default : None
        A pseudo random number generator used for the initialization of the
        residuals when eigen_solver == 'arpack'.

    Attributes
    ----------

    lambdas_ :
        Eigenvalues of the centered kernel matrix

    alphas_ :
        Eigenvectors of the centered kernel matrix

    dual_coef_ :
        Inverse transform matrix

    X_transformed_fit_ :
        Projection of the fitted data on the kernel entropy components

    References
    ----------
    Kernel ECA based on:
    (c) Robert Jenssen, University of Tromso, Norway, 2010 
        R. Jenssen, "Kernel Entropy Component Analysis,"
        IEEE Trans. Patt. Anal. Mach. Intel., 32(5), 847-860, 2010.

    """

    def __init__(self, n_components=None, kernel="linear",
                 gamma=None, degree=3, coef0=1, kernel_params=None, eigen_solver='auto',
                 tol=0, max_iter=None, random_state=None,center=False):
        self.n_components = n_components
        self._kernel = kernel
        self.kernel_params = kernel_params
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.eigen_solver = eigen_solver #not selectable right now
        self.tol = tol
        self.max_iter = max_iter
        self.random_state = random_state
        self._centerer = KernelCenterer()
        self.center = center
        
    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def _get_kernel(self, X, Y=None):
        if callable(self._kernel):
            params = self.kernel_params or {}
        else:
            params = {"gamma": self.gamma,
                      "degree": self.degree,
                      "coef0": self.coef0}
        return pairwise_kernels(X, Y, metric=self._kernel,
                                filter_params=True, **params)

    def _fit_transform(self, K):
        """ Fit's using kernel K"""
        # center kernel
        if self.center == True:
            K = self._centerer.fit_transform(K)

        X_transformed = self.kernelECA(K=K)    
        self.X_transformed = X_transformed
        return K        

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
        K = self._get_kernel(X)
        self._fit_transform(K)

        self.X_fit_ = X
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

        X_transformed= self.X_transformed
        
        return X_transformed

    def transform(self, X):
        """Transform X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """
        check_is_fitted(self, 'X_fit_')

        K = self._centerer.transform(self._get_kernel(X, self.X_fit_))
        return np.dot(K, self.alphas_ / np.sqrt(self.lambdas_))

    def inverse_transform(self, X):
        raise NotImplementedError("Function inverse_transform is not implemented.")

    # here are the helper functions => to integrate in the code!
    def kernelECA(self,K):   
        if self.n_components is None:
            n_components = K.shape[0]
        else:
            n_components = min(K.shape[0], self.n_components)
             
        # compute eigenvectors
        self.lambdas_, self.alphas_ = linalg.eigh(K)
        
        d = self.lambdas_
        E = self.alphas_
        # sort eigenvectors in descending order
        D,E = self.sort_eigenvalues(d,E)
    
        d = np.diag(D)
        sorted_entropy_index,entropy = self.ECA(D,E)
        Es = E[:,sorted_entropy_index]
        ds = d[sorted_entropy_index]

        Phi = np.zeros((K.shape[0],n_components))
        for i in range(n_components):
            Phi[:,i] = np.sqrt(ds[i]) * Es[:,i]
    
        X_transformed = Phi
        
        return X_transformed

    def sort_eigenvalues(self,D,E):
        d = D
        indices = np.argsort(d)[::-1]
    
        d = d[indices]
        D = np.zeros((len(d),len(d)))
        for i in range(len(d)):
            D[i,i] = d[i]
        E = E[:,indices]

        return D,E
    
    def ECA(self,D,E):
        N = E.shape[0]
        entropy = np.multiply(np.diag(D).T , (np.dot(np.ones((1,N)),E))**2)[0]
        indices = np.argsort(entropy)[::-1]
        entropy = entropy[indices]
        return indices,entropy