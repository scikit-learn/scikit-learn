import numpy as np
import scipy.linalg

from scikits.learn.utils.utils import fast_logdet

def bayesian_ridge( X , Y, step_th=300,th_w = 1.e-12,ll_bool=False) :
    """
    Bayesian ridge regression. Optimize the regularization parameter alpha
    within a simple bayesian framework (MAP).

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
    w : numpy array of shape (dim)
         mean of the weights distribution.
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
    See Bishop p 345-348 for more details.
    """

    beta = 1./np.var(Y)
    alpha = 1.0

    log_likelihood = []
    has_converged = False
    gram = np.dot(X.T, X)
    ones = np.eye(gram.shape[1])
    sigma = scipy.linalg.pinv(alpha*ones + beta*gram)
    w = np.dot(beta*sigma,np.dot(X.T,Y))
    old_w = w
    while not has_converged and step_th:

	### Update Parameters
	# alpha
        lmbd_ = np.real(scipy.linalg.eigvals(beta * gram.T))
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
	if ll_bool :
	  residual_ = (Y - np.dot(X, w))**2
	  ll = 0.5*X.shape[1]*np.log(alpha) + 0.5*X.shape[0]*np.log(beta)
	  ll -= (0.5*beta*residual_.sum()+ 0.5*alpha*np.dot(w.T,w))
	  ll -= fast_logdet(alpha*ones + beta*gram)
	  ll -= X.shape[0]*np.log(2*np.pi)
	  log_likelihood.append(ll)

    return w, log_likelihood


def bayesian_linear(alpha, beta):
    """
    Like bayesian_ridge,
    but alpha, beta is given
    """

    ### Compute mu and sigma
    gram = np.dot(X.T, X)
    ones = np.eye(gram.shape[1])
    sigma = scipy.linalg.pinv(alpha*ones + beta*gram)
    w = np.dot(beta*sigma,np.dot(X.T,Y))


    return w



class BayesianRegression(object):
    """
    Encapsulate various bayesian regression algorithms
    """
    
    def __init__(self):
        pass

    def fit(self, X, Y):
        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)
        self.w,self.log_likelihood = bayesian_ridge(X, Y)

    def predict(self, T):
        return np.dot(T, self.w)


