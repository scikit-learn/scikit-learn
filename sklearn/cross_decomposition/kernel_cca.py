"""Kernel Canonical Correlation Analysis (KCCA)"""

# Author: Huibin Shen <icdishb@gmail.com>
# License: BSD 3 clause

import cmath
import warnings
import numpy as np
from scipy import linalg

from ..utils import check_array
from ..utils.arpack import eigsh
from ..base import BaseEstimator
from sklearn.preprocessing import KernelCenterer
from ..metrics.pairwise import pairwise_kernels

__all__ = ['KernelCCA']

class KernelCCA(BaseEstimator):
    """Kernel Canonical Correlation Analysis (KCCA)

    Canonical correlation analysis between two kernels. By default, 
    the kernel matrices are decomposed by Partial Gram-Schmidt 
    Orthgonalization. Kernels are centered before everything.

    Parameters
    ----------
    n_components: int, optional
        Number of components (default is 2),

    kernel: {'linear', 'poly', 'rbf', 'sigmoid', 
             'cosine', 'precomputed'}, optional
        Kernel types (default is 'linear').

    degree : int, optional
        Degree for poly kernels (default is 3).
        Ignored by other kernels.

    gamma : float, optional
        Kernel coefficient for rbf and poly kernels (Default: 1/n_features).
        Ignored by other kernels.

    coef0 : float, optional
        Independent term in poly and sigmoid kernels.
        Ignored by other kernels.

    kernel_params : mapping of string to any, optional
        Parameters (keyword arguments) and values for kernel passed as
        callable object. Ignored by other kernels.

    eigen_solver : {'auto', 'dense', 'arpack'}, optional
        Select eigensolver to use.  If n_components is much less than
        the number of training samples, arpack may be more efficient
        than the dense eigensolver (default is 'auto').

    center : boolean, optional
        Center the data ? If the data is already centered, set center to
        False to not center it internally again (default is False).

    pgso : boolean, optional
        Use partial Gram-Schmidt orthogonalization to decompose kernel
        matrix (default is True).

    eta : float, optional
        Precision parameter used in partial Gram-Schmidt 
        orthogonalization (default is 0).

    kapa : float, optional
        Regulation parameter in kernel CCA, 
        values between 0 to 1 (default is 0.1).

    nor : {1, 2, 3}, optional
        Normalize option (default is 2).
        nor=1 without normalisation and use Gram-Schmidt space.
        nor=2 normalisation in kernel space.
        nor=3 normalisation in Gram-Schmidt space.

    tol: float, optional
        Convergence tolerance for arpack (default is 0, 
        optimal value will be chosen by arpack).

    max_iter : int, optional
        Maximum number of iterations for arpack (default is None,
        optimal value will be chosen by arpack)

    copy : boolean, optional
        Whether a forced copy will be triggered (default is False, 
        used in check_arrays).


    Attributes
    ----------
    KXc_ : array_like, shape = (n_samples, n_components)
        Centered kernel matrix for X or KX. If your data is already 
        centered, KXc_ is just the kernel matrix without centering again.

    KYc_ : array_like, shape = (n_samples, n_components)
        Centered kernel matrix for Y or KY. If your data is already 
        centered, KYc_ is just the kernel matrix without centering again.

    Rx_ : array_like, shape = (n_samples, k)
        If pgso is true, the kernel matrix is approximated by KXc_ = 
        Rx_ * RX_'. If pgso is false, then Rx_ is None.

    Ry_ : array_like, shape = (n_samples, k)
        If pgso is true, the kernel matrix is approximated by KYc_ = 
        Ry_ * Ry_'. If pgso is false, then Ry_ is None.

    lambdas_ : array_like, shape = (n_components,)
        Canonical correlation coefficients

    alphas_ : array_like, shape = (n_samples, n_components)
        Basis matrix (components) for KXc.

    betas_ : array_like, shape = (n_samples, n_components)
        Basis matrix (components) for KYc.


    Notes
    -----
    Find a set of basis vectors alphas (for KX) and betas (for KY) 
    such that the corr(KX*alphas[:,i], KY*betas[:,i]) is maximized,
    where i = 0 : n_compoenents. Kernels are centered before everything.

    By default, Partial Gram-Schmidt Orthogonalization (pgso) is used to 
    decompose the kernel matrix. When PGSO is not used, only positive 
    correlation coefficients and corresponding basis vectors  are returned.

    This implementation is based on the matlab code for kernel CCA provided
    by David Hardoon. 
    http://www.davidroihardoon.com/Professional/Code_files/kcca_package.tar.gz

    Examples
    --------
    >>> from sklearn.cross_decomposition import KernelCCA
    >>> X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [3.,5.,4.]]
    >>> Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]
    >>> kcca = KernelCCA(kernel='linear', n_components=2, kapa=0.1)
    >>> kcca.fit(X, Y)
    ... # doctest: +NORMALIZE_WHITESPACE
    KernelCCA(center=False, coef0=1, copy=True, degree=3, eigen_solver='auto', 
         eta=0, gamma=None, kapa=0.1, kernel='linear', kernel_params=None,     
         max_iter=500, n_components=2, nor=2, pgso=True, tol=1e-06) 
    >>> kcca.lambdas_
    array([ 0.9997788 ,  0.76983395])
    >>> kcca.alphas_
    array([[-0.08682032,  1.36279428],
           [ 0.01811873,  0.13090947],
           [-0.00897666,  0.50482967],
           [ 0.02814122, -0.36253905]])
    >>> kcca.betas_
    array([[-0.02072478,  3.25968106],
           [ 0.01283091, -1.98755628],
           [-0.02694432,  4.43065658],
           [ 0.01536388, -2.03628425]])

    References
    ----------

    David R. Hardoon, Sandor Szedmak, John Shawe-Taylor. 2004.
    Canonical Correlation Analysis: An Overview with Application to
    Learning Models. Neural computation, volume 16, no. 12,
    pp2639-2664, 2004. 


    See also
    --------
    CCA
    """

    def __init__(self, n_components=2, kernel="linear", gamma=None, 
                 degree=3, coef0=1, kernel_params=None, eigen_solver='auto',
                 center=False, pgso=True, eta=0, kapa=0.1, nor=2, 
                 max_iter=500, tol=1e-6, copy=True):
        self.n_components = n_components
        self.kernel = kernel
        self.kernel_params = kernel_params
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.eigen_solver = eigen_solver
        self.center = center
        self.pgso = pgso
        self.eta = eta
        self.kapa = kapa
        self.nor = nor
        self.max_iter = max_iter
        self.tol= tol
        self.copy = copy


    def _get_kernel(self, X, Y=None):
        if callable(self.kernel):
            params = self.kernel_params or {}
        else:
            params = {"gamma": self.gamma,
                      "degree": self.degree,
                      "coef0": self.coef0}
        return pairwise_kernels(X, Y, metric=self.kernel,
                                filter_params=True, **params)


    def _solve_eigenvalues(self, A):
        """Solve eigenvalue problem for matrix A"""

        n = A.shape[0]

        if self.eigen_solver == 'auto':
            if n > 200 and self.n_components < 10:
                eigen_solver = 'arpack'
            else:
                eigen_solver = 'dense'
        else:
            eigen_solver = self.eigen_solver

        if eigen_solver == 'dense':
            rr, alphas = linalg.eigh(
                A, eigvals=(n - self.n_components, n - 1))
        elif eigen_solver == 'arpack':
            rr, alphas = eigsh(A, self.n_components, which="LA", 
                               tol=self.tol,maxiter=self.max_iter)
        return (rr, alphas)


    def _pgso(self, K, eta):
        """Partial Gram Schmidt Orthogonalization"""

        m = K.shape[0]
        index = np.zeros((m,1))
        tSize = np.zeros((m,1))
        feat = np.zeros((m,m))
        
        norm2 = np.diag(K)
        j = 0
        while np.sum(norm2) > eta and j != m:
            # find the best new element
            j2 = np.argmax(norm2)
            
            # save the index and setting
            index[j,0] = j2
            tSize[j,0] = np.sqrt(norm2[j2])

            # calculate new features
            tSum = np.dot(feat[0:m,0:j], feat[j2,0:j].T)
            feat[0:m,j] = (K[0:m,j2] - tSum)/tSize[j]

            # updating diagonal elements
            norm2 = norm2 - feat[0:m,j]**2
            j += 1
        feat = feat[0:m,0:j]
        return (feat, tSize, index)


    def _normalisekcca(self, A, K):
        """Normalise features"""

        n = A.shape[0]
        KK = np.dot(K.T, K)

        comp = np.sqrt(np.diag(np.dot(np.dot(A.T, KK), A)))
        comp[comp==0] = 1
        comp = comp.reshape((A.shape[1],1))

        re = A / np.dot(np.ones((n,1)), comp.T)
        return re


    def _fit(self, Kx, Ky):
        """Kernel CCA without PGSO."""

        I = np.diag(np.ones(Kx.shape[0]))
        Kxinv = linalg.pinv(np.dot((1-self.kapa)*I,Kx) + self.kapa*I)
        Kyinv = linalg.pinv(np.dot((1-self.kapa)*I,Ky) + self.kapa*I)
        A = np.dot(np.dot(np.dot(Kxinv, Ky), Kyinv), Kx)
        A = 0.5*(A.T + A) + np.diag(np.ones(A.shape[0]))*10e-6

        rankA = np.linalg.matrix_rank(A)
        if self.n_components > rankA:
            warnings.warn("set n_components to %d" % rankA)
            self.n_components = rankA
        rr, alpha = self._solve_eigenvalues(A) 

        #if np.all(rr.real > 0): # rr is lambda^2 
        #    r = np.sqrt(rr.real)
        #else: # rr has negative eigenvalues
        #    r = np.array(map(cmath.sqrt, rr.real))

        # throw away negative eigenvalues
        if not np.all(rr.real > 0):
            alpha = alpha[:,rr>0]
            rr = rr[rr>0]
            self.n_components = len(rr)
            warnings.warn("set n_components to %d" % len(rr))
        r = np.sqrt(rr.real)

        # recover lambda
        beta = np.dot(np.dot(Kyinv, Kx), alpha)
        t = Ky.shape[0]
        beta = beta / np.tile(r, (t,1))
                    
        # correct correlation, alphas and betas
        nalpha = self._normalisekcca(alpha, Kx)
        nbeta = self._normalisekcca(beta, Ky)
        r = np.diag(np.dot(np.dot(np.dot(nalpha.T, Kx.T), Ky), nbeta))
        
        # sort eigenvalues in descending order
        inds = np.argsort(r)[::-1]
        self.lambdas_ = r[inds]
        self.alphas_ = nalpha[:,inds]
        self.betas_ = nbeta[:,inds]
        self.Rx_ = None
        self.Ry_ = None


    def _fit_pgso(self, Kx, Ky):
        """Kernel CCA with PGSO on the kernel matrices."""

        # Decompose kernel with PGSO
        Rx, Rxsize, Rxindex = self._pgso(Kx, self.eta)
        Ry, Rysize, Ryindex = self._pgso(Ky, self.eta)
        self.Rx_ = Rx
        self.Ry_ = Ry

        # Creating covariance matrix with the new feature matrix Rx, Ry
        Zxx = np.dot(Rx.T, Rx)
        Zxy = np.dot(Rx.T, Ry)
        Zyy = np.dot(Ry.T, Ry)
        Zyx = Zxy.T
        tEyeY = np.diag(np.ones(Zyy.shape[0]))
        tEyeX = np.diag(np.ones(Zxx.shape[0]))

        # solve the eigen value problem
        B = np.dot((1-self.kapa)*tEyeX, Zxx) + self.kapa*tEyeX
        S = linalg.cholesky(B, lower=True)
        invS = linalg.pinv(S)
        Zyyinv = linalg.pinv(np.dot((1-self.kapa)*tEyeY,Zyy) + self.kapa*tEyeY)
        SinvZxy = np.dot(invS, Zxy)
        A = np.dot(np.dot(SinvZxy, Zyyinv), SinvZxy.T)
        A = 0.5*(A.T + A) + np.diag(np.ones(A.shape[0]))*10e-6

        rankA = np.linalg.matrix_rank(A)
        if self.n_components > rankA:
            warnings.warn("set n_components to %d" % rankA)
            self.n_components = rankA
        rr, pha = self._solve_eigenvalues(A) # pha is \hat{alpha} in the paper

        # recover lambda
        if np.all(rr.real > 0): # rr is lambda^2
            r = np.sqrt(rr.real)
        else:
            r = np.array(map(cmath.sqrt, rr.real))

        # recover alpha and beta
        alpha = np.dot(invS.T, pha) # alpha is \tilde{alpha} in the paper
        invRx = np.dot(Rx, linalg.pinv(np.dot(Rx.T, Rx)))
        nalpha = np.dot(invRx, alpha) # nalpha is alpha in the paper

        beta = np.dot(np.dot(Zyyinv, Zyx), alpha) # beta is \tilde{beta}
        t = Zyy.shape[0]
        beta = beta / np.tile(r, (t,1))
        invRy = np.dot(Ry, linalg.pinv(np.dot(Ry.T, Ry)))
        nbeta = np.dot(invRy, beta) # nbeta is beta in the paper

        # normalising the feature vectors and recomputing the correlation
        if self.nor == 1: # no normalizsation but using GS projected space
            nalpha = alpha
            nbeta = beta
            warnings.warn("Remember to project your testing Kernels into the GS space!")
        elif self.nor == 2: # normalisation in kernel space
            nalpha = self._normalisekcca(nalpha, Kx)
            nbeta = self._normalisekcca(nbeta, Ky)
            r = np.diag(np.dot(np.dot(np.dot(nalpha.T, Kx.T), Ky), nbeta))
        elif self.nor == 3: # normalisation in GS space
            nalpha = self._normalisekcca(alpha, Rx)
            nbeta = self._normalisekcca(beta, Ry)
            r = np.diag(np.dot(np.dot(np.dot(nalpha.T, Rx.T), Ry), nbeta))
            warnings.warn("Remember to project your testing Kernels into the GS space!")
        else:
            warnings.warn("nor can only be 1, 2 or 3")
        
        # sort eigenvalues in descending order
        inds = np.argsort(r)[::-1]
        self.lambdas_ = r[inds]
        self.alphas_ = nalpha[:,inds]
        self.betas_ = nbeta[:,inds]


    def fit(self, X, Y):
        """Fit the KCCA model with two views represented by kernels X and Y.

        Parameters
        ----------
        X : array_like, shape = (n_samples, n_features) for data matrix
            or shape = (n_samples, n_samples) for kernel matrix.
            When both X and Y are kernel matrix, the kernel parameter 
            should be set to 'precomputed'.
            It is considered to be one view of the data.

        Y : array_like, shape = (n_samples, n_features) for data matrix
            or shape = (n_samples, n_samples) for kernel matrix.
            When both X and Y are kernel matrix, the kernel parameter 
            should be set to 'precomputed'.
            It is considered to be another view of the data.
        
        Returns
        -------
        self : object
            Returns the instance itself.
        """
        X = check_array(X, dtype=np.float, copy=self.copy)
        Y = check_array(Y, dtype=np.float, copy=self.copy)

        n = X.shape[0]
        p = X.shape[1]
        q = Y.shape[1]

        if n != Y.shape[0]:
            raise ValueError(
                'Incompatible shapes: X has %s samples, while Y '
                'has %s' % (X.shape[0], Y.shape[0]))
        if self.n_components < 1 or self.n_components > n:
            raise ValueError('Invalid number of components')
        if self.eigen_solver not in ("auto", "dense", "arpack"):
            raise ValueError("Got eigen_solver %s when only 'auto', "
                             "'dense' and 'arparck' are valid" % 
                             self.algorithm)
        if self.kernel == 'precomputed' and (p != n or q != n):
            raise ValueError('Invalid kernel matrices dimension')
        if not self.pgso and (self.kapa <= 0 or self.kapa >= 1):
            raise ValueError('kapa should be in (0, 1) when pgso=False')
        if self.pgso and (self.kapa < 0 or self.kapa > 1):
            raise ValueError('kapa should be in [0, 1] when pgso=True')


        KX = self._get_kernel(X)
        KY = self._get_kernel(Y)
        
        if self.center:
            kc = KernelCenterer()
            self.KXc_ = kc.fit_transform(KX)
            self.KYc_ = kc.fit_transform(KY)
        else:
            self.KXc_ = KX
            self.KYc_ = KY

        if self.pgso: # use PGSO to decompose kernel matrix
            self._fit_pgso(self.KXc_, self.KYc_)
        else:
            self._fit(self.KXc_, self.KYc_)
        return self
