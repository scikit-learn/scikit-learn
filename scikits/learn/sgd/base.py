# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Stochastic Gradient Descent (SGD) with sparse data. """

import numpy as np
from scipy import sparse

from ..base import BaseEstimator


class LinearModel(BaseEstimator):
    """Linear Model trained by minimizing a regularized training
    error using SGD.
    
    Parameters
    ----------
    loss : Loss
        The loss function to be used. 
    penalty : str, ('l2'|'l1'|'elasticnet')
        The penalty (aka regularization term) to be used.
    alpha : float
        Constant that multiplies the penalty. Defaults to 0.0001.
    rho : float
        The ElasticNet mixing parameter, with 0 < rho <= 1. 
    coef_ : ndarray of shape n_features
        The initial coeffients to warm-start the optimization
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

    def __init__(self, penalty='l2', loss="hinge", alpha = 0.0001,
		 rho = 0.85, coef_ = None, intercept_ = 0.0,
		 fit_intercept = True):
        self.penalty = penalty
        self.loss = loss
	self.alpha = alpha
        self.rho = rho
	self.coef_ = coef_
	self.intercept_ = intercept_
        self.fit_intercept = fit_intercept
	self._get_penalty_type()

    def _get_penalty_type(self):
	penalty_types = {"l2":2, "l1":1, "elasticnet":3}
	try:
	    self.penalty_type = penalty_types[self.penalty]
	except KeyError:
	    raise ValueError("The penalty %s is not supported. " % self.penalty)

    def predict(self, X):
	"""Predict using the linear model

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Returns predicted class labeles (either 1 or -1 or 0).
        """
	X = np.asanyarray(X)
	# FIXME what if prediction is 0 - break randomly?
        return np.sign(np.dot(X, self.coef_) + self.intercept_)

    def predict_margin(self, T):
	"""Predict signed 'distance' to the hyperplane (aka confidence score). 

	If 

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Returns signed 'distance' to the hyperplane
        """
	X = np.asanyarray(X)
        return np.dot(X, self.coef_) + self.intercept_

    def predict_proba(self, T):
        # how can this be, logisitic *does* implement this
        raise NotImplementedError(
                'sgd does not provide this functionality')
