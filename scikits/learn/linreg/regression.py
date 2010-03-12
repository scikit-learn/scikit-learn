# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

# $Id$
"""
Regression using linear models: Linear, Ridge.
"""

import numpy as np
import scipy.linalg

class LinearRegression(object):
    """
    Linear Regression.

    Parameters
    ----------
    This class takes no parameters

    Members
    -------
    coef_ : array
        Estimated coefficients for the linear regression problem.

    This is just plain linear regression wrapped is a Predictor object.
    """

    def fit(self,X,Y):
        """
        Fit linear model
        """
        self.coef_, self.residues_, self.rank_, self.singular_ = \
                scipy.linalg.lstsq(X, Y)
        return self

    def predict(self, T):
        """
        Predict using linear model
        """
        return np.dot(T, self.coef_)


class Ridge(object):
    """
    Ridge regression.


    Parameters
    ----------
    alpha : ridge parameter. Small positive values of alpha improve
    the coditioning of the problem and reduce the variance of the
    estimates.
    
    Examples
    --------
    # With more samples than features
    >>> import numpy as np
    >>> nsamples, nfeatures = 10, 5
    >>> np.random.seed(0)
    >>> Y = np.random.randn(nsamples)
    >>> X = np.random.randn(nsamples, nfeatures)
    >>> clf = Ridge(alpha=alpha)
    >>> clf.fit(X, Y)
    ?

    See also
    --------
    http://scikit-learn.sourceforge.net/doc/modules/linreg.html
    """

    def __init__(self, alpha=1.0):
        self.alpha = alpha

    def fit(self, X, Y):
        """Fit Ridge regression model"""
        nsamples, nfeatures = X.shape

        if nsamples > nfeatures:
            # w = inv(X^t X + alpha*Id) * X.T y
            self.coef_ = scipy.linalg.solve(np.dot(X.T,X) + self.alpha * np.eye(nfeatures),
                                  np.dot(X.T, Y))
        else:
            # w = X.T * inv(X X^t + alpha*Id) y
            self.coef_ = np.dot(X.T,
                    scipy.linalg.solve(np.dot(X, X.T) + self.alpha * np.eye(nsamples), Y))

        return self

    def predict(self, T):
        """
        Predict using Linear Model
        """
        return np.dot(T, self.coef_)


class BayesianRidge:
    """
    Encapsulate various bayesian regression algorithms
    """

    def __init__(self, ll_bool=False, step_th=300, th_w=1.e-12):
        self.ll_bool = ll_bool
        self.step_th = step_th
        self.th_w = th_w

    def fit(self, X, Y):
        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)
        self.w, self.alpha, self.beta, self.sigma, self.log_likelihood = \
            bayesian_regression_ridge(X, Y, self.step_th, self.th_w, self.ll_bool)
        return self

    def predict(self, T):
        return np.dot(T, self.w)


### helper methods
### we should homogeneize this

def bayesian_regression_ridge( X , Y, step_th=300, th_w = 1.e-12, ll_bool=False):
    """
    Bayesian ridge regression. Optimize the regularization parameters alpha
    (precision of the weights) and beta (precision of the noise) within a simple
    bayesian framework (MAP).

    Parameters
    ----------
    X : numpy array of shape (length,features)
    data
    Y : numpy array of shape (length)
    target
    step_th : int (defaut is 300)
          Stop the algorithm after a given number of steps.
    th_w : float (defaut is 1.e-12)
       Stop the algorithm if w has converged.
    ll_bool  : boolean (default is False).
           If True, compute the log-likelihood at each step of the model.

    Returns
    -------
    w : numpy array of shape (nb_features)
         mean of the weights distribution.
    alpha : float
       precision of the weights.
    beta : float
       precision of the noise.
    sigma : numpy array of shape (nb_features,nb_features)
        variance-covariance matrix of the weights
    log_likelihood : list of float of size steps.
             Compute (if asked) the log-likelihood of the model.

    Examples
    --------
    >>> X = np.array([[1], [2]])
    >>> Y = np.array([1, 2])
    >>> w = ridge_regression(X,Y)
    w = 1.

    Notes
    -----
    See Bishop p 167-169 for more details.
    """

    beta = 1./np.var(Y)
    alpha = 1.0

    log_likelihood = []
    has_converged = False
    gram = np.dot(X.T, X)
    ones = np.eye(gram.shape[1])
    sigma = scipy.linalg.pinv(alpha*ones + beta*gram)
    w = np.dot(beta*sigma,np.dot(X.T,Y))
    old_w = np.copy(w)
    eigen = np.real(scipy.linalg.eigvals(gram.T))
    while not has_converged and step_th:

        ### Update Parameters
        # alpha
        lmbd_ = np.dot(beta, eigen)
        gamma_ = (lmbd_/(alpha + lmbd_)).sum()
        alpha = gamma_/np.dot(w.T, w)

        # beta
        residual_ = (Y - np.dot(X, w))**2
        beta = (X.shape[0]-gamma_) / residual_.sum()

        ### Compute mu and sigma
        sigma = scipy.linalg.pinv(alpha*ones + beta*gram)
        w = np.dot(beta*sigma,np.dot(X.T,Y))
        step_th -= 1

        # convergence : compare w
        has_converged =  (np.sum(np.abs(w-old_w))<th_w)
        old_w = w

    ### Compute the log likelihood
    if ll_bool:
        residual_ = (Y - np.dot(X, w))**2
        ll = 0.5*X.shape[1]*np.log(alpha) + 0.5*X.shape[0]*np.log(beta)
        ll -= (0.5*beta*residual_.sum()+ 0.5*alpha*np.dot(w.T,w))
        ll -= fast_logdet(alpha*ones + beta*gram)
        ll -= X.shape[0]*np.log(2*np.pi)
        log_likelihood.append(ll)

    return w, alpha, beta, sigma, log_likelihood
