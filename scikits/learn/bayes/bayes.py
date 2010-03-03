import numpy as np
import scipy.linalg


def fast_logdet(A):
    """
    Compute log(det(A)) for A symmetric
    Equivalent to : np.log(nl.det(A))
    but more robust
    It returns -Inf if det(A) is non positive
    Copyright : A. Gramfort 2010   
    """
    from math import exp,log
    ld = np.sum(np.log(np.diag(A)))
    if not np.isfinite(ld):
        return -np.inf
    a = exp(ld/A.shape[0])
    d = scipy.linalg.det(A/a)
    if d <= 0:
        return -np.inf
    ld += log(d)
    if not np.isfinite(ld):
        return -np.inf
    return ld



def bayesian_ridge( X , Y, step_th=300,th_w = 1.e-6,ll_bool=True) :
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
    th_w : float (defaut is 1.e-6)
	   Stop the algorithm if w has converged.
    ll_bool  : boolean (default is True).
	       If True, compute the log-likelihood at each step of the model.

    Returns
    -------
    w : numpy array of shape (dim)
         mean of the weights distribution.
   
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
	  ll = 0.5*X.shape[1]*np.log(alpha) + 0.5*X.shape[0]*np.log(beta)
	  ll -= 0.5*beta*residual_.sum()+ 0.5*alpha*np.dot(w.T,w)
	  ll -= fast_logdet(inv_sigma_) 
	  ll -= X.shape[0]*np.log(2*np.pi)
	  log_likelihood.append(ll)

    return w

