""" Partial Least Square
"""

# Author: Edouard Duchesnay <edouard.duchesnay@cea.fr>
# License: BSD Style.

#from .base import BaseEstimator
from scikits.learn.base import BaseEstimator
import warnings
import numpy as np
from scipy import linalg

def _nipals_twoblocks_inner_loop(X, Y, mode="A", max_iter = 500, tol = 1e-06):
    """Inner loop of the iterative NIPALS algorithm. provide an alternative
    of the svd(X'Y) ie. return the first left and rigth singular vectors of X'Y
    See PLS for the meaning of the parameters.
    """
    p = X.shape[1]
    q = Y.shape[1]
    u = np.ones((p,1))/p
    v = np.ones((q,1))/q
    
    u_old = 0
    ite = 1
    # Inner loop of the Wold algo.
    while True:
        print " ite",ite
        # Update y_score: the Y latent scores
        y_score = np.dot(Y, v)/np.dot(v.T, v)
        
        # Update u: the X loadings
        u = np.dot(X.T, y_score)/np.dot(y_score.T, y_score)
        u = u / np.sqrt(np.dot(u.T,u))
        
        # Update x_score: the X latent scores
        x_score    = np.dot(X, u)/np.dot(u.T, u)
        
        # Update v: the Y loadings
        v = np.dot(Y.T, x_score)/np.dot(x_score.T, x_score)
        v = v / np.sqrt(np.dot(v.T,v))
        
        u_diff = u - u_old
        if np.dot(u_diff.T, u_diff) < tol:
            break
        if ite == max_iter:
            warnings.warn('Maximum number of iterations reached')
            break
        u_old = u
        #v_old = v
        ite += 1
    return u, v

def _svd_cross_product(X, Y):
    C = np.dot(X.T,Y)
    U, s, Vh = linalg.svd(C, full_matrices = False)
    u = U[:,[0]]
    v = Vh.T[:,[0]]
    return u, v
    
class PLS(BaseEstimator):
    """Partial Least Square (PLS)

    We use the therminology defined by [Wegelin et al. 2000]
    This implementation uses the PLS-W2A (PLS Wold 2 blocks) algorithm as
    referenced by Wegelin et al. 2000 which is based on two nested loop:
    The outer loop iterate over compements. 
    The inner loop is a two block generalization of the of power method 
    usually use to compute eigenvectors on singlese block.

    Parameters; ,k 
    ----------
    X: array-like of predictors, shape (n_samples, p)
        Training vector, where n_samples in the number of samples and
        p is the number of predictors.

    Y: array-like of response, shape (n_samples, q)
        Training vector, where n_samples in the number of samples and
        q is the number of response variables.

    n_components: int, number of components to keep. (default 2).

    deflation_mode: str, "canonical" or "regression". "canonical" performs
        a symetric deflation (each block is deflated on its own latent
        variable: X on x_score and Y on y_score). "regression" performs an asymetric
        deflation: both X and Y are deflated on x_score. The two modes yield to same
        results on the first component, differences start with the second
        components. The canonical mode is mostly used for modeling, while
        regression is prediction oriented.
        
    mode: "A" classical PLS and "B" CCA.
    
    center_X: boolean, center X? (default True)
    
    center_Y: boolean, center Y? (default True)
        
    scale_X: boolean, scale X? (default True)
            
    scale_X: boolean, scale X? (default True)
    
    algorithm: str "nipals" or "svd" the algorithm used to estimate the 
        loadings, it will be called "n_components" time ie.: for each iteration
        of the outer loop.
    
    max_iter: an integer, the maximum number of iterations (default 500) of the
        NIPALS inner loop (used only if algorithm="nipals")
    
    tol: a not negative real, the tolerance used in the iterative algorithm
         default 1e-06.
 
    Notes
    -----
    PLS mode A, regression, ie.: the PLS regression also known as PLS2
    
    PLS mode A, canonical a symetric version of the PLS regression
    
        Optimize:
        max corr(Xu, Yv) * var(Xu), |u| = |v| = 1
    

    
    PLS mode C, canonical, ie.: the CCA
    
    References
    ----------
    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with 
    emphasis on the two-block case. Technical Report 371, Department of 
    Statistics, University of Washington, Seattle, 2000.
 
    Attributes
    ----------
    x_loadings_: array, [p x n_components] loadings for the X block
    y_loadings_: array, [q x n_components] loadings for the Y block
    x_score_: array, [p x n_samples] scores for the X block
    y_score_: array, [q x n_samples] scores for the Y block


    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn.pls import PLS
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> Y = np.array([[-3, -2], [-1, -2], [-5, -3], [4, 1], [-2, 2], [5, 6]])
    >>> pls = PLS(n_components=2)
    >>> pls.fit(X, Y)
    PCA(copy=True, n_components=2, whiten=False)
    >>> print pca.explained_variance_ratio_
    [ 0.99244289  0.00755711]

    See also
    --------
    PLS_SVD
    CCA

    """
    def __init__(self, n_components=2, deflation_mode = "canonical", mode = "A",
                 center_X = True, center_Y = True,
                 scale_X  = True, scale_Y  = True,
                 algorithm = "nipals",
                 max_iter = 500, tol = 1e-06):
        self.n_components = n_components
        self.deflation_mode = deflation_mode
        self.mode = mode
        self.center_X = center_X
        self.center_Y = center_Y
        self.scale_X  = scale_X
        self.scale_Y  = scale_Y
        self.algorithm = algorithm
        self.max_iter = max_iter
        self.tol      = tol

    def fit(self, X, Y,  **params):
        self._set_params(**params)
        X = np.asanyarray(X)
        y = np.asanyarray(Y)

        n = X.shape[0]
        p = X.shape[1]
        q = Y.shape[1]

        if X.ndim != 2:
            raise ValueError('X must be a 2D array')

        if n != Y.shape[0]:
            raise ValueError(
                'Incompatible shapes: X has %s samples, while Y '
                'has %s' % (X.shape[0], Y.shape[0]))

        if self.n_components < 1 or self.n_components > p:
            raise ValueError('invalid number of components')
            
        if self.algorithm is "svd" and self.mode is "B":
            raise ValueError(
                'Incompatible configuration: mode B is not implemented with svd'
                'algorithm')
        if not self.deflation_mode in ["canonical","regression"]:
            raise ValueError(
                'The deflation mode is not valid')

        # Scale the data
        if self.center_X:
            self.X_mu_  = X.mean(axis=0)
            X = (X - self.X_mu_)  
        if self.scale_X :
            self.X_std_ = X.std(axis=0,ddof=1)
            X = X / self.X_std_
        if self.center_Y:
            self.Y_mu_  = Y.mean(axis=0)
            Y = Y - self.Y_mu_  
        if self.scale_Y :
            self.Y_std_ = Y.std(axis=0,ddof=1)
            Y = Y / self.Y_std_
        # Residuals (deflated) matrices
        Xk = X
        Yk = Y
        # Results matrices
        self.x_scores_ = np.zeros((n, self.n_components))
        self.y_scores_ = np.zeros((n, self.n_components))
        self.x_loadings_ = np.zeros((p, self.n_components))
        self.y_loadings_ = np.zeros((q, self.n_components))
        self.x_regs_     = np.zeros((p, self.n_components))#Xk'x_scorek/x_scorek'x_scorek
        # NIPALS algo: outer loop, over components
        for k in xrange(self.n_components):
            #print " comp",k,"--------------------------------------------------"
            #1) Loadings estimation (inner loop)
            # -----------------------------------
            if self.algorithm is "nipals":
                u, v = _nipals_twoblocks_inner_loop(
                        X=Xk, Y=Yk, mode=self.mode, 
                        max_iter=self.max_iter, tol=self.tol)
            if self.algorithm is "svd":
                u, v = _svd_cross_product(X=Xk, Y=Yk)
            # compute scores
            x_score = np.dot(Xk, u)
            y_score = np.dot(Yk, v)
            #2) Deflation
            # ------------
            # - regress Xk's on x_score
            a = np.dot(Xk.T, x_score)/np.dot(x_score.T, x_score) # p x 1
            # - substract rank-one approximations to obtain remainder matrix
            Xk = Xk - np.dot(x_score, a.T)
            if self.deflation_mode is "canonical":
                # - regress Yk's on y_score, then substract rank-one approximation
                b  = np.dot(Yk.T, y_score)/np.dot(y_score.T, y_score) # q x 1
                Yk = Yk - np.dot(y_score, b.T)
            if self.deflation_mode is "regression":
                # - regress Yk's on x_score, then substract rank-one approximation
                b  = np.dot(Yk.T, x_score)/np.dot(x_score.T, x_score) # q x 1
                Yk = Yk - np.dot(x_score, b.T)
            # 3) Store loadings and scores
            self.x_scores_[:,k] = x_score.ravel()
            self.y_scores_[:,k] = y_score.ravel()
            self.x_loadings_[:,k] = u.ravel()
            self.y_loadings_[:,k] = v.ravel()
            self.x_regs_[:,k]     = a.ravel()
            
            print "component",k,"----------------------------------------------"
            print "X rank-one approximations and residual"
            print np.dot(x_score, a.T)
            print Xk
            print "Y rank-one approximations and residual"
            print np.dot(y_score, b.T)
            print Yk
        return self



class PLS_SVD(BaseEstimator):
    """Partial Least Square SVD
    
    Simply perform a svd on the crosscovariance matrix: X'Y
    The are no iterative deflation here.

    Parameters
    ----------
    X: array-like of predictors, shape (n_samples, p)
        Training vector, where n_samples in the number of samples and
        p is the number of predictors.

    Y: array-like of response, shape (n_samples, q)
        Training vector, where n_samples in the number of samples and
        q is the number of response variables.

    n_components: int, number of components to keep. (default 2).
    
    center_X: boolean, center X? (default True)
    
    center_Y: boolean, center Y? (default True)
        
    scale_X: boolean, scale X? (default True)
            
    scale_X: boolean, scale X? (default True)
     
    Attributes
    ----------
    u_: array, [p x n_components] loadings for the X block
    v_: array, [q x n_components] loadings for the Y block
    x_score_: array, [p x n_samples] scores for X the block
    y_score_: array, [q x n_samples] scores for the Y block

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn.pls import PLS

    See also
    --------
    PLS
    CCA

    """
    def __init__(self, n_components=2,
                 center_X = True, center_Y = True,
                 scale_X  = True, scale_Y  = True):
        self.n_components = n_components
        self.center_X = center_X
        self.center_Y = center_Y
        self.scale_X  = scale_X
        self.scale_Y  = scale_Y
 
    def fit(self, X, Y,  **params):
        self._set_params(**params)
        X = np.asanyarray(X)
        y = np.asanyarray(Y)
        
        n = X.shape[0]
        p = X.shape[1]
        q = Y.shape[1]

        if X.ndim != 2:
            raise ValueError('X must be a 2D array')

        if n != Y.shape[0]:
            raise ValueError(
                'Incompatible shapes: X has %s samples, while Y '
                'has %s' % (X.shape[0], Y.shape[0]))

        if self.n_components < 1 or self.n_components > p:
            raise ValueError('invalid number of components')

        # Scale the data
        if self.center_X:
            self.X_mu_  = X.mean(axis=0)
            X -= self.X_mu_  
        if self.scale_X :
            self.X_std_ = X.std(axis=0)
            X /= self.X_std_
        if self.center_Y:
            self.Y_mu_  = Y.mean(axis=0)
            Y -= self.Y_mu_  
        if self.scale_Y :
            self.Y_std_ = Y.std(axis=0)
            Y /= self.Y_std_
        
        C = np.dot(X.T,Y)
        U, s, V = linalg.svd(C, full_matrices = False)
        V = V.T
        self.x_score_    = np.dot(X, U)
        self.y_score_ = np.dot(Y, V)
        self.u_     = U
        self.v_     = V
        return self


