#from .base import BaseEstimator
from scikits.learn.base import BaseEstimator
import warnings
import numpy as np

class PLS(BaseEstimator):
    def __init__(self, n_components=2, mode = "canonical", ):
        self.n_components = n_components
        self.mode = mode

    def fit(self, X, Y, max_iter = 500, tol = 1e-06,
            center_X = True, center_Y = True, scale_X = True, scale_Y = True,
            **params):
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
        if center_X:
            self.X_mu_  = X.mean(axis=0)
            X -= self.X_mu_  
        if scale_X :
            self.X_std_ = X.std(axis=0)
            X /= self.X_std_
        if center_Y:
            self.Y_mu_  = Y.mean(axis=0)
            Y -= self.Y_mu_  
        if scale_Y :
            self.Y_std_ = Y.std(axis=0)
            Y /= self.Y_std_
        # Residuals (deflated) matrices
        X_r = X
        Y_r = Y
        # Results matrices
        self.xi_    = np.zeros((n, self.n_components))
        self.omega_ = np.zeros((n, self.n_components))
        self.u_     = np.zeros((p, self.n_components))
        self.v_     = np.zeros((q, self.n_components))

        # Outer loop: over components
        for r in xrange(self.n_components):
            print " comp",r
            u = np.ones((p,1))/p
            v = np.ones((q,1))/q
            u_old = 0
            v_old = 0
            k = 1
            # Inner loop of the Wold algo.
            while True:
                print " ite",k
                #k += 1
                # Update the latent scores
                xi    = np.dot(X_r, u)
                omega = np.dot(Y_r, v)
                
                # Update loadings
                u = np.dot(X_r.T, omega)/np.dot(omega.T, omega)
                u = u / np.sqrt(np.dot(u.T,u))
                
                v = np.dot(Y_r.T, xi)/np.dot(xi.T, xi)
                v = v / np.sqrt(np.dot(v.T,v))
                
                u_diff = u - u_old
                
                if np.dot(u_diff.T, u_diff) < tol:
                    break
                if k == max_iter:
                    warnings.warn('Maximum number of iterations reached for component %s' % (h))
                    break
                u_old = u
                #v_old = v
                k = k+1
            # Deflation: substract rank-one approximations to obtanin remainder matrices
            X_r = X_r - np.dot(xi, np.dot(xi.T, X_r) / np.dot(xi.T, xi))
            if self.mode is "canonical":
                Y_r = Y_r - np.dot(omega, np.dot(omega.T, Y_r) / np.dot(omega.T, omega))   
            if self.mode is "regression":
                Y_r = Y_r - np.dot(xi, np.dot(xi.T, Y_r) / np.dot(xi.T, xi))
                
            self.xi_[:,r]    = xi.ravel()
            self.omega_[:,r] = omega.ravel()
            self.u_[:,r]     = u.ravel()
            self.v_[:,r]     = v.ravel()
        return self
