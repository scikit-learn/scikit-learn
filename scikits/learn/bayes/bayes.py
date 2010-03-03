import numpy as np
import scipy.linalg

def bayesian_ridge( X , Y, step_th=300,th_w = 1.e-6, verbose = True) :
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
    verbose  : boolean (default is True).
	       Set the output on the console.

    Returns
    -------
    w : numpy array of shape (dim)
         mean of the weights distribution.
   


    Notes
    -----
    See Bishop p 345-348 for more details.
    """

    beta = 1./np.var(Y)
    alpha = 1.0

    has_converged = False
    gram = np.dot(X.T, X)
    ones = np.eye(gram.shape[1])
    sigma = scipy.linalg.pinv(alpha*ones + beta*gram)
    w = np.dot(beta*sigma,np.dot(X.T,Y))
    while not has_converged and step_th:

	w_ = np.copy(w)
       
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
	if np.sum(np.abs(w-w_))<th_w :
	    break

    return w

