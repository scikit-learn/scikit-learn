import numpy as np
import scipy.linalg
from scikits.learn.utils.utils import fast_logdet


def bayesian_regression_noprior( X , Y, ll_bool=False):
    """
    Bayesian regression. Find the solutions by Maximum-likelihood (ML) -
    equivalent to the OLS estimate with a noise precision beta.

    Parameters
    ----------
    X : numpy array of shape (length,features)
	data
    Y : numpy array of shape (length)
	target
    ll_bool  : boolean (default is False).
	       If True, compute the log-likelihood at each step of the model.

    Returns
    -------
    w : numpy array of shape (nb_features)
         mean of the weights distribution.
    beta : float
	   precision of the noise.
    log_likelihood : list of float of size steps.
		     Compute (if asked) the log-likelihood of the model.
   
    Examples
    --------

    Notes
    -----
    See Bishop p 138-143 for more details.
    """
   
    gram = np.dot(X.T, X)
    w = np.dot(scipy.linalg.pinv(gram),np.dot(X.T,Y))
    residual_ = (Y - np.dot(X, w))**2
    beta = X.shape[0] / residual_.sum()
    log_likelihood = []
    if ll_bool :
      log_likelihood = 0.5*X.shape[0]*np.log(beta)
      log_likelihood -= 0.5*X.shape[0]*np.log(2*np.pi) 
      log_likelihood -= 0.5*beta*residual_.sum()
    return w,beta,log_likelihood






def bayesian_regression_ridge( X , Y, step_th=300,th_w = 1.e-12,ll_bool=False) :
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

    return w,alpha,beta,sigma,log_likelihood



def bayesian_regression_ard(X , Y, step_th=300, th_w=1.e-12,\
			    alpha_th=1.e+16, ll_bool=False):
    """
    Bayesian ard-based regression. Optimize the regularization parameters alpha
    (vector of precisions of the weights) and beta (precision of the noise).


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
    alpha_th : number
           threshold on the alpha, to avoid divergence. Remove those features
	   from the weights computation if is alpha > alpha_th  (default is
	    1.e+16).
    ll_bool  : boolean (default is False).
	       If True, compute the log-likelihood at each step of the model.

    Returns
    -------
    w : numpy array of shape (nb_features)
         mean of the weights distribution.
    alpha : numpy array of shape (nb_features)
	   precision of the weights.
    beta : float
	   precision of the noise.
    sigma : numpy array of shape (nb_features,nb_features)
	    variance-covariance matrix of the weights
    log_likelihood : list of float of size steps.
		     Compute (if asked) the log-likelihood of the model.
   
    Examples
    --------

    Notes
    -----
    See Bishop p 345-348 for more details.
    """
    gram = np.dot(X.T, X)
    beta = 1./np.var(Y)
    alpha = np.ones(gram.shape[1])

    
    log_likelihood = None
    if ll_bool :
      log_likelihood = []
    has_converged = False
    ones = np.eye(gram.shape[1])
    sigma = scipy.linalg.pinv(alpha*ones + beta*gram)
    w = np.dot(beta*sigma,np.dot(X.T,Y))
    old_w = np.copy(w)
    keep_a  = np.ones(X.shape[1],dtype=bool)
    while not has_converged and step_th:

	
	# alpha
	gamma_ = 1 - alpha[keep_a]*np.diag(sigma)
	alpha[keep_a] = gamma_/w[keep_a]**2
  
	# beta
	residual_ = (Y - np.dot(X[:,keep_a], w[keep_a]))**2
	beta = (X.shape[0]-gamma_.sum()) / residual_.sum()

	### Avoid divergence of the values by setting a maximum values of the
	### alpha
	keep_a = alpha<alpha_th
	gram = np.dot(X.T[keep_a,:], X[:,keep_a])

        ### Compute mu and sigma
	ones = np.eye(gram.shape[1])
	sigma = scipy.linalg.pinv(alpha[keep_a]*ones+ beta*gram)
	w[keep_a] = np.dot(beta*sigma,np.dot(X.T[keep_a,:],Y))
        step_th -= 1

	# convergence : compare w
	has_converged =  (np.sum(np.abs(w-old_w))<th_w)
        old_w = w

	
	### Compute the log likelihood
	if ll_bool :
	  A_ = np.eye(X.shape[1])/alpha
	  C_ = (1./beta)*np.eye(X.shape[0]) + np.dot(X,np.dot(A_,X.T))
	  ll = X.shape[0]*np.log(2*np.pi)+fast_logdet(C_)
	  ll += np.dot(Y.T,np.dot(scipy.linalg.pinv(C_),Y)) 
	  log_likelihood.append(-0.5*ll)
	
    return w,alpha,beta,sigma,log_likelihood




class BayesianRegression(object):
    """
    Encapsulate various bayesian regression algorithms
    """
    def __init__(self, ll_bool=False):
	self.ll_bool = ll_bool

    def fit(self,X,Y):
        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)
	self.w, self.beta, self.log_likelihood = \
	bayesian_regression_noprior(X, Y, self.ll_bool)
  
    def predict(self, T):
        return np.dot(T, self.w)


class RidgeRegression(BayesianRegression):
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
	self.w ,self.alpha ,self.beta ,self.sigma ,self.log_likelihood = \
	bayesian_regression_ridge(X, Y, self.step_th, self.th_w, self.ll_bool)



class ARDRegression(BayesianRegression):
    """
    Encapsulate various bayesian regression algorithms
    """
    
    def __init__(self, ll_bool=False, step_th=300, th_w=1.e-12,\
		alpha_th=1.e+16):
        self.ll_bool = ll_bool
	self.step_th = step_th
	self.th_w = th_w
	self.alpha_th = alpha_th

    def fit(self, X, Y):
        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)
	self.w ,self.alpha ,self.beta ,self.sigma ,self.log_likelihood = \
	bayesian_regression_ard(X, Y, self.step_th, self.th_w,\
			    self.alpha_th, self.ll_bool)
