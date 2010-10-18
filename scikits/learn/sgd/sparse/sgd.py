# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Implementation of Stochastic Gradient Descent (SGD) with sparse data."""

import warnings
import numpy as np
from scipy import sparse

from ...base import ClassifierMixin
from ..base import LinearModel
from . import sgd_fast_sparse

class SGD(LinearModel, ClassifierMixin):
    """Linear Model trained by minimizing a regularized training
    error using SGD.

    This implementation works on scipy.sparse X and dense coef_.
    
    Parameters
    ----------
    loss : Loss
        The loss function to be used. 
    penalty : str, ('l2'|'l1'|'elasticnet')
        The penalty (aka regularization term) to be used.
    reg : float
        The regularization parameter, i.e. tradeoff between loss and penalty.
    alpha : float
        Constant that multiplies the L1 term. Defaults to 1.0
    rho : float
        The ElasticNet mixing parameter, with 0 < rho <= 1.
    coef_ : ndarray of shape n_features
        The initial coeffients to warm-start the optimization
    intercept_ : float
        The initial intercept to warm-start the optimization
    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        Weights asigned to the features.

    `intercept_` : float, shape = [n_class-1]
        Constants in decision function.

    """
    
    def _set_coef(self, coef_):
        self.coef_ = coef_
        if coef_ is None:
            self.sparse_coef_ = None
        else:
            n_features = len(coef_)
            # sparse representation of the fitted coef for the predict method
            self.sparse_coef_ = sparse.csr_matrix(coef_)

    def fit(self, X, Y, n_iter=5, **params):
        """Fit current model with SGD

        X is expected to be a sparse matrix. For maximum efficiency, use a
        sparse matrix in CSR format (scipy.sparse.csr_matrix)
        """
        self._set_params(**params)
        X = sparse.csr_matrix(X)
        Y = np.asanyarray(Y, dtype = np.float64)

        n_samples, n_features = X.shape[0], X.shape[1]
        if self.coef_ is None:
            self.coef_ = np.zeros(n_features, dtype=np.float64, order = "c")
        
        X_data = np.array(X.data, np.float32, order = "c")
	X_indices = X.indices
	X_indptr = X.indptr
	loss = sgd_fast_sparse.Hinge()
	verbose = 2
	shuffle = 0
	coef_, intercept_ = sgd_fast_sparse.plain_sgd(self.coef_, self.intercept_,
						      loss, int(self.penalty_type),
						      self.alpha, self.rho, X_data,
						      X_indices, X_indptr, Y,
						      int(n_iter),
						      int(self.fit_intercept),
						      int(verbose), int(shuffle))

        # update self.coef_ and self.sparse_coef_ consistently
        self._set_coef(coef_)
	self.intercept_ = intercept_

        # return self for chaining fit and predict calls
        return self

    ## def predict(self, X):
## 	"""Predict using the linear model

##         Parameters
##         ----------
##         X : numpy array of shape [n_samples, n_features]

##         Returns
##         -------
##         C : array, shape = [n_samples]
##             Returns predicted class labeles (either 1 or -1 or 0).
##         """
## 	raise NotImplemented("not implemented yet")

##     def predict_margin(self, T):
## 	"""Predict signed 'distance' to the hyperplane (aka confidence score). 

## 	If 

##         Parameters
##         ----------
##         X : numpy array of shape [n_samples, n_features]

##         Returns
##         -------
##         C : array, shape = [n_samples]
##             Returns signed 'distance' to the hyperplane
##         """
## 	raise NotImplemented("not implemented yet")
