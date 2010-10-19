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
    loss : str, ('hinge'|'log'|'modifiedhuber')
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

    def _get_loss_function(self):
	loss_functions = {"hinge" : sgd_fast_sparse.Hinge(),
			  "log" : sgd_fast_sparse.Log(),
			  "modifiedhuber" : sgd_fast_sparse.ModifiedHuber(),
			  }
	try:
	    self.loss_function = loss_functions[self.loss]
	except KeyError:
	    raise ValueError("The loss %s is not supported. " % self.loss)
    
    def _set_coef(self, coef_):
        self.coef_ = coef_
        if coef_ is None:
            self.sparse_coef_ = None
        else:
            n_features = len(coef_)
            # sparse representation of the fitted coef for the predict method
            self.sparse_coef_ = sparse.csr_matrix(coef_)

    def fit(self, X, Y, n_iter = 5, **params):
        """Fit current model with SGD

        X is expected to be a sparse matrix. For maximum efficiency, use a
        sparse matrix in CSR format (scipy.sparse.csr_matrix)
        """
        self._set_params(**params)
	self.n_iter = int(n_iter)
        X = sparse.csr_matrix(X)
        Y = np.asanyarray(Y, dtype = np.float64)

        n_samples, n_features = X.shape[0], X.shape[1]
        if self.coef_ is None:
            self.coef_ = np.zeros(n_features, dtype=np.float64, order = "c")
        
        X_data = np.array(X.data, np.float64, order = "c")
	X_indices = X.indices
	X_indptr = X.indptr
	verbose = 2
	shuffle = 0
	print "norm: ", self.penalty_type
	print "rho: ", self.rho
	print "alpha:", self.alpha
	coef_, intercept_ = sgd_fast_sparse.plain_sgd(self.coef_, self.intercept_,
						      self.loss_function, self.penalty_type,
						      self.alpha, self.rho, X_data,
						      X_indices, X_indptr, Y,
						      self.n_iter, int(self.fit_intercept),
						      verbose, int(shuffle))

        # update self.coef_ and self.sparse_coef_ consistently
        self._set_coef(coef_)
	self.intercept_ = intercept_

        # return self for chaining fit and predict calls
        return self

    def predict(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] with the predicted class labels (either -1, 1, or 0).
        """        
        return np.sign(self.predict_margin(X))

    def predict_margin(self, X):
        """Predict signed 'distance' to the hyperplane (aka confidence score). 

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] with signed 'distances' to the hyperplane.
        """
        # np.dot only works correctly if both arguments are sparse matrices
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        return np.ravel(np.dot(self.sparse_coef_, X.T).todense()
                        + self.intercept_)
