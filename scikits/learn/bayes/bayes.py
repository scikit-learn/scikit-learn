import numpy as np
import scipy.linalg

def bayesian_regression( X , Y, step_th=300, verbose = True) :
    """
    Regularized regression function.
    Function for fitting an Automatic Relevance Determination to the data.
        See Bishop p 345-348 for more details.

    Parameters
    ----------
    phi : numpy array of shape (length,dim)
          functionnal images (gram matrix)
    y : numpy array of shape (length)
        target.
    prune_th : number
           Defaut is 1.e+12. If not None, we remove the alpha by
           removing the ones supra-thresholded.
    mu_th : number
        threshold on the delta of the weights to stop the convergence.
        Defaut is 1.e-4.
    lambd_th : number
           threshold on the lambda, to avoid divergence. Set the lambda
           to lambda_th (default is 1.e+12) is lambda > lambda_th.
    step_th : number.
          Stop the algorithm if the number of step is > step_th.
    mode : string
           mode of computing for alpha : direct differenciation
           "DirectDiff" (defaut)(see Bishop p347), or
           expectation-maximization "EM" (see Bishop p450).
           "DirectDiff" is normally faster.
    verbose  : boolean.
           Set the output on the console (default is True).

    Returns
    -------
    mu : numpy array of shape (dim)
         mean of the weights distribution.
    log_evidence : number
               the log evidence of p(y|0,Sigma)    
    """

    beta = 1./np.var(Y)
    alpha = 1.0

    has_converged = False
    gram = np.dot(X.T, X)
    ones = np.eye(gram.shape[1])
    sigma = scipy.linalg.pinv(alpha*ones + beta*gram)
    w = np.dot(beta*sigma,np.dot(X.T,Y))
    while not has_converged and step_th:

        ### Update Parameters
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

    return w

