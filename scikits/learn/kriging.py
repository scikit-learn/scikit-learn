#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    This module contains a contribution to the scikit-learn project that
    implements Gaussian Process based prediction (also known as Kriging).
"""
__authors__ = u"Vincent Dubourg <vincent.dubourg@gmail.com>"
__date__ =  u"10/31/2010"

################
# Dependencies #
################

import numpy as np
from numpy import random
from scipy import optimize
from numpy import linalg
from .base import BaseEstimator

#########################
# Convenience functions #
#########################

def colsum(x):
    """
    Performs columnwise sums of elements in x depending on its shape.
    
    Parameters
    ----------
    x : array_like
	An array of size (p, q).
    
    Returns
    -------
    s : array_like
	An array of size (q, ) which contains the columnwise sums of the elements in x.
    """
    
    x = np.asanyarray(x, dtype=np.float)
    
    if x.ndim > 1:
	s = x.sum(axis=0)
    else:
	s = x
    
    return s


#################################
# Defaut autocorrelation models #
#################################

def correxp1(theta,d):
    """
    Exponential autocorrelation model.
    (Ornstein-Uhlenbeck stochastic process)
						   n
    correxp1 : theta, dx --> r(theta, dx) = exp(  sum  - theta_i * |dx_i| )
					 	 i = 1
    Parameters
    ----------
    theta : array_like
	An array with shape 1 (isotropic) or n (anisotropic) giving the autocorrelation parameter(s).
    dx : array_like
	An array with shape (n_eval, n_features) giving the componentwise distances between locations x and x' at which the correlation model should be evaluated.
    
    Returns
    -------
    r : array_like
	An array with shape (n_eval, ) containing the values of the autocorrelation model.
    """
    
    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)
    
    if d.ndim > 1:
	n_features = d.shape[1]
    else:
	n_features = 1
    if theta.size == 1:
	theta = np.repeat(theta, n_features)
    elif theta.size != n_features:
	raise ValueError, "Length of theta must be 1 or "+str(n_features)
    
    td = - theta.reshape(1, n_features) * abs(d)
    r = np.exp(np.sum(td,1))
    
    return r

def correxp2(theta,d):
    """
    Squared exponential correlation model (Radial Basis Function).
    (Infinitely differentiable stochastic process, very smooth)
						   n
    correxp2 : theta, dx --> r(theta, dx) = exp(  sum  - theta_i * (dx_i)^2 )
						 i = 1
    Parameters
    ----------
    theta : array_like
	An array with shape 1 (isotropic) or n (anisotropic) giving the autocorrelation parameter(s).
    dx : array_like
	An array with shape (n_eval, n_features) giving the componentwise distances between locations x and x' at which the correlation model should be evaluated.
    
    Returns
    -------
    r : array_like
	An array with shape (n_eval, ) containing the values of the autocorrelation model.
    """
    
    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)
    
    if d.ndim > 1:
	n_features = d.shape[1]
    else:
	n_features = 1
    if theta.size == 1:
	theta = np.repeat(theta, n_features)
    elif theta.size != n_features:
	raise Exception, "Length of theta must be 1 or "+str(n_features)
  
    td = - theta.reshape(1, n_features) * d**2
    r = np.exp(np.sum(td,1))
    
    return r

def correxpg(theta,d):
    """
    Generalized exponential correlation model.
    (Useful when one does not know the smoothness of the function to be predicted.)
						   n
    correxpg : theta, dx --> r(theta, dx) = exp(  sum  - theta_i * |dx_i|^p )
						 i = 1
    Parameters
    ----------
    theta : array_like
	An array with shape 1+1 (isotropic) or n+1 (anisotropic) giving the autocorrelation parameter(s) (theta, p).
    dx : array_like
	An array with shape (n_eval, n_features) giving the componentwise distances between locations x and x' at which the correlation model should be evaluated.
    
    Returns
    -------
    r : array_like
	An array with shape (n_eval, ) with the values of the autocorrelation model.
    """
    
    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)
    
    if d.ndim > 1:
	n_features = d.shape[1]
    else:
	n_features = 1
    lth = theta.size
    if n_features > 1 and lth == 2:
	theta = np.hstack([np.repeat(theta[0], n_features), theta[1]])
    elif lth != n_features+1:
	raise Exception, "Length of theta must be 2 or "+str(n_features+1)
    else:
	theta = theta.reshape(1, lth)
    
    td = - theta[:,0:-1].reshape(1, n_features) * abs(d)**theta[:,-1]
    r = np.exp(np.sum(td,1))
    
    return r

def corriid(theta,d):
    """
    Spatial independence correlation model (pure nugget).
    (Useful when one wants to solve an ordinary least squares problem!)
						    n
    corriid : theta, dx --> r(theta, dx) = 1 if   sum |dx_i| == 0 and 0 otherwise
						  i = 1
    Parameters
    ----------
    theta : array_like
	None.
    dx : array_like
	An array with shape (n_eval, n_features) giving the componentwise distances between locations x and x' at which the correlation model should be evaluated.
    
    Returns
    -------
    r : array_like
	An array with shape (n_eval, ) with the values of the autocorrelation model.
    """
    
    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)
    
    n_eval = d.shape[0]
    r = np.zeros(n_eval)
    # The ones on the diagonal of the correlation matrix are enforced within
    # the KrigingModel instanciation to allow multiple design sites in this
    # ordinary least squares context.
    
    return r

def corrcubic(theta,d):
    """
    Cubic correlation model.
					       n
    corrcubic : theta, dx --> r(theta, dx) = prod max(0, 1 - 3(theta_j*d_ij)^2 + 2(theta_j*d_ij)^3) ,  i = 1,...,m
					     j = 1
    Parameters
    ----------
    theta : array_like
	An array with shape 1 (isotropic) or n (anisotropic) giving the autocorrelation parameter(s).
    dx : array_like
	An array with shape (n_eval, n_features) giving the componentwise distances between locations x and x' at which the correlation model should be evaluated.
    
    Returns
    -------
    r : array_like
	An array with shape (n_eval, ) with the values of the autocorrelation model.
    """
    
    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)
    
    if d.ndim > 1:
	n_features = d.shape[1]
    else:
	n_features = 1
    lth = theta.size
    if  lth == 1:
	theta = np.repeat(theta, n_features)[np.newaxis,:]
    elif lth != n_features:
	raise Exception, "Length of theta must be 1 or "+str(n_features)
    else:
	theta = theta.reshape(1, n_features)
    
    td = abs(d) * theta
    td[td > 1.] = 1.
    ss = 1. - td**2. * (3. - 2.*td)
    r = np.prod(ss,1)
    
    return r

def corrlin(theta,d):
    """
    Linear correlation model.
					     n
    corrlin : theta, dx --> r(theta, dx) = prod max(0, 1 - theta_j*d_ij) ,  i = 1,...,m
					   j = 1
    Parameters
    ----------
    theta : array_like
	An array with shape 1 (isotropic) or n (anisotropic) giving the autocorrelation parameter(s).
    dx : array_like
	An array with shape (n_eval, n_features) giving the componentwise distances between locations x and x' at which the correlation model should be evaluated.
    
    Returns
    -------
    r : array_like
	An array with shape (n_eval, ) with the values of the autocorrelation model.
    """
    
    theta = np.asanyarray(theta, dtype=np.float)
    d = np.asanyarray(d, dtype=np.float)
    
    if d.ndim > 1:
	n_features = d.shape[1]
    else:
	n_features = 1
    lth = theta.size
    if  lth == 1:
	theta = np.repeat(theta, n_features)[np.newaxis,:]
    elif lth != n_features:
	raise Exception, "Length of theta must be 1 or "+str(n_features)
    else:
	theta = theta.reshape(1, n_features)
    
    td = abs(d) * theta
    td[td > 1.] = 1.
    ss = 1. - td
    r = np.prod(ss,1)
	
    return r


############################
# Defaut regression models #
############################

def regpoly0(x):
    """
    Zero order polynomial (constant, p = 1) regression model.
    
    regpoly0 : x --> f(x) = 1
    
    Parameters
    ----------
    x : array_like
	An array with shape (n_eval, n_features) giving the locations x at which the regression model should be evaluated.
    
    Returns
    -------
    f : array_like
	An array with shape (n_eval, p) with the values of the regression model.
    """
    
    x = np.asanyarray(x, dtype=np.float)
    n_eval = x.shape[0]
    f = np.ones([n_eval,1])
    
    return f

def regpoly1(x):
    """
    First order polynomial (hyperplane, p = n) regression model.
    
    regpoly1 : x --> f(x) = [ x_1, ..., x_n ].T
    
    Parameters
    ----------
    x : array_like
	An array with shape (n_eval, n_features) giving the locations x at which the regression model should be evaluated.
    
    Returns
    -------
    f : array_like
	An array with shape (n_eval, p) with the values of the regression model.
    """
    
    x = np.asanyarray(x, dtype=np.float)
    n_eval = x.shape[0]
    f = np.hstack([np.ones([n_eval,1]), x])
    
    return f

def regpoly2(x):
    """
    Second order polynomial (hyperparaboloid, p = n*(n-1)/2) regression model.
    
    regpoly2 : x --> f(x) = [ x_i*x_j,  (i,j) = 1,...,n ].T
					    i > j
    
    Parameters
    ----------
    x : array_like
	An array with shape (n_eval, n_features) giving the locations x at which the regression model should be evaluated.
    
    Returns
    -------
    f : array_like
	An array with shape (n_eval, p) with the values of the regression model.
    """
    
    x = np.asanyarray(x, dtype=np.float)
    n_eval, n_features = x.shape
    f = np.hstack([np.ones([n_eval,1]), x])
    for  k in range(n_features):
	f = np.hstack([f, x[:,k,np.newaxis] * x[:,k:]])
    
    return f


###########################
# The Kriging Model class #
###########################

class KrigingModel(BaseEstimator):
    """
    A class that implements scalar Gaussian Process based prediction (also known as Kriging).

    Example
    -------
    import numpy as np
    from scikits.learn.kriging import KrigingModel, correxp2
    import pylab as pl

    f = lambda x: x*np.sin(x)
    X = np.array([1., 3., 5., 6., 7., 8.])
    Y = f(X)
    aKrigingModel = KrigingModel(corr=correxp2)
    aKrigingModel.fit(X, Y, theta0=1e-1, thetaL=1e-3, thetaU=1e0, random_start=100)

    pl.figure(1)
    x = np.linspace(0,10,1000)
    y, MSE = aKrigingModel.predict(x, eval_MSE=True)
    sigma = np.sqrt(MSE)
    pl.plot(x, f(x), 'r:', label=u'$f(x) = x\,\sin(x)$')
    pl.plot(X, Y, 'r.', markersize=10, label=u'Observations')
    pl.plot(x, y, 'k-', label=u'$\widehat{f}(x) = {\\rm BLUP}(x)$')
    pl.fill(np.concatenate([x, x[::-1]]), np.concatenate([y - 1.9600 * sigma, (y + 1.9600 * sigma)[::-1]]), alpha=.5, fc='b', ec='None', label=u'95\% confidence interval')
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.legend(loc='upper left')

    pl.figure(2)
    theta_values = np.logspace(np.log10(aKrigingModel.thetaL),np.log10(aKrigingModel.thetaU),100)
    score_values = []
    for t in theta_values:
	    score_values.append(aKrigingModel.score(theta=t)[0])
    pl.plot(theta_values, score_values)
    pl.xlabel(u'$\\theta$')
    pl.ylabel(u'Score')
    pl.xscale('log')

    pl.show()

    Methods
    -------
    fit(X, y) : self
        Fit the model.

    predict(X) : array
        Predict using the model.

    Todo
    ----
    o Add the 'sparse' storage mode for which the correlation matrix is stored in its sparse eigen decomposition format instead of its full Cholesky decomposition.

    Notes
    -----
    See references:
    [1] H.B. Nielsen, S.N. Lophaven, H. B. Nielsen and J. Sondergaard (2002). DACE - A MATLAB Kriging Toolbox.
	http://www2.imm.dtu.dk/~hbn/dace/dace.pdf
    """
    
    def __init__(self, regr = regpoly0, corr = correxp2, beta0 = None, storage_mode = 'full', verbose = True):
	"""
	The Kriging Model constructor.

	Parameters
	----------
	regr : lambda function, optional
	    A regression function returning an array of outputs of the linear regression functional basis.
	    (The number of observations m should be greater than the size p of this basis)
	    Default assumes a simple constant regression trend (see regpoly0).
	    
	corr : lambda function, optional
	    A stationary autocorrelation function returning the autocorrelation between two points x and x'.
	    Default assumes a squared-exponential autocorrelation model (see correxp2).
	    
	beta0 : double array_like, optional
	    Default assumes Universal Kriging (UK) so that the vector beta of regression weights is estimated
	    by Maximum Likelihood. Specifying beta0 overrides estimation and performs Ordinary Kriging (OK).
	    
	storage_mode : string, optional
	    A string specifying whether the Cholesky decomposition of the correlation matrix should be
	    stored in the class (storage_mode = 'full') or not (storage_mode = 'light').
	    Default assumes storage_mode = 'full', so that the Cholesky decomposition of the correlation matrix is stored.
	    This might be a useful parameter when one is not interested in the MSE and only plan to estimate the BLUP,
	    for which the correlation matrix is not required.
	    
	verbose : boolean, optional
	    A boolean specifying the verbose level.
	    Default is verbose = True.

	Returns
	-------
	aKrigingModel : self
	    A Kriging Model object awaiting data to be fitted to.
	"""
	
	self.regr = regr
	self.corr = corr
	self.beta0 = beta0
	
	# Check storage mode
	if storage_mode != 'full' and storage_mode != 'light':
	    if storage_mode == 'sparse':
		raise NotImplementedError, "The 'sparse' storage mode is not supported yet. Please contribute!"
	    else:
		raise ValueError, "Storage mode should either be 'full' or 'light'. Unknown storage mode: "+str(storage_mode)
	else:
	    self.storage_mode = storage_mode
	
	self.verbose = verbose
	
	return None
    
    def fit(self, X, y, theta0 = 1e-1, thetaL = None, thetaU = None, optimizer = 'fmin_cobyla', random_start = 1, normalize = True, nugget = 1e-15):
	"""
	The Kriging Model fitting method.

	Parameters
	----------
	X : double array_like
	    An array with shape (n_samples, n_features) with the design sites at which observations were made.
	    
	y : double array_like
	    An array with shape (n_features, ) with the observations of the scalar output to be predicted.
	    
	theta0 : double array_like, optional
	    An array with shape (n_features, ) or (1, ).
	    The parameters in the autocorrelation model.
	    If thetaL and thetaU are also specified, theta0 is considered as the starting point
	    for the Maximum Likelihood Estimation of the best set of parameters.
	    Default assumes isotropic autocorrelation model with theta0 = 1e-1.
	    
	thetaL : double array_like, optional
	    An array with shape matching theta0's.
	    Lower bound on the autocorrelation parameters for Maximum Likelihood Estimation.
	    Default is None, skips Maximum Likelihood Estimation and uses theta0.
	    
	thetaU : double array_like, optional
	    An array with shape matching theta0's.
	    Upper bound on the autocorrelation parameters for Maximum Likelihood Estimation.
	    Default is None and skips Maximum Likelihood Estimation and uses theta0.
	    
	normalize : boolean, optional
	    Design sites X and observations y are centered and reduced wrt means and standard deviations
	    estimated from the n_samples observations provided.
	    Default is true so that data is normalized to ease MLE.
	    
	nugget : double, optional
	    Introduce a nugget effect to allow smooth predictions from noisy data.
	    Default assumes no nugget effect (nugget = 0.).
	    
	optimizer : string, optional
	    A string specifying the optimization algorithm to be used ('fmin_cobyla' is the sole algorithm implemented yet).
	    Default uses 'fmin_cobyla' algorithm from scipy.optimize.
	    
	random_start : int, optional
	    The number of times the Maximum Likelihood Estimation should be performed from a random starting point.
	    The first MLE always uses the specified starting point (theta0),
	    the next starting points are picked at random according to an exponential distribution (log-uniform in [thetaL, thetaU]).
	    Default does not use random starting point (random_start = 1).
	
	Returns
	-------
	aKrigingModel : self
	    A fitted Kriging Model object awaiting data to perform predictions.
	"""
	
	# Force data to numpy.array type (from coding guidelines)
        X = np.asanyarray(X, dtype=np.float)
        y = np.asanyarray(y, dtype=np.float)
        
	# Check design sites & observations
	n_samples_X = X.shape[0]
	if X.ndim > 1:
	    n_features = X.shape[1]
	else:
	    n_features = 1
	X = X.reshape(n_samples_X, n_features)
	
	n_samples_y = y.shape[0]
	if y.ndim > 1:
	    raise NotImplementedError, "y has more than one dimension. This is not supported yet (scalar output prediction only). Please contribute!"
	y = y.reshape(n_samples_y, 1)
	
	if n_samples_X != n_samples_y:
	    raise Exception, "X and y must have the same number of rows!"
	else:
	    n_samples = n_samples_X
	
	# Check correlation parameters
	theta0 = np.asanyarray(theta0, dtype=np.float)
	if theta0.ndim == 0:
	    theta0 = theta0[np.newaxis]
	lth = theta0.size
	if thetaL != None and thetaU != None: # Parameters optimization case
	    thetaL = np.asanyarray(thetaL, dtype=np.float)
	    if thetaL.ndim == 0:
		thetaL = thetaL[np.newaxis]
	    thetaU = np.asanyarray(thetaU, dtype=np.float)
	    if thetaU.ndim == 0:
		thetaU = thetaU[np.newaxis]
	    if thetaL.size != lth or thetaU.size != lth:
		raise ValueError, "theta0, thetaL and thetaU must have the same length"
	    if np.any(thetaL <= 0) or np.any(thetaU < thetaL):
		raise ValueError, "The bounds must satisfy O < thetaL <= thetaU"
	elif thetaL == None and thetaU == None: # Given parameters case
	    if np.any(theta0 <= 0):
		raise ValueError, "theta0 must be strictly positive"
	elif thetaL == None or thetaU == None: # Error
	    raise Exception, "thetaL and thetaU should either be both or neither specified"
	
	# Normalize data or don't
	if normalize:
	    mean_X = np.mean(X, axis=0)
	    std_X = np.sqrt(1./(n_samples-1.)*np.sum((X - np.mean(X,0))**2.,0)) #np.std(y, axis=0)
	    mean_y = np.mean(y, axis=0)
	    std_y = np.sqrt(1./(n_samples-1.)*np.sum((y - np.mean(y,0))**2.,0)) #np.std(y, axis=0)
	    std_X[std_X == 0.] = 1.
	    std_y[std_y == 0.] = 1.
	else:
	    mean_X = np.array([0.])
	    std_X = np.array([1.])
	    mean_y = np.array([0.])
	    std_y = np.array([1.])
	
	X_ = (X - mean_X) / std_X
	y_ = (y - mean_y) / std_y
	
	# Calculate matrix of distances D between samples
	mzmax = n_samples*(n_samples-1)/2
	ij = np.zeros([mzmax, 2])
	D = np.zeros([mzmax, n_features])
	ll = np.array([-1])
	for k in range(n_samples-1):
	    ll = ll[-1] + 1 + range(n_samples-k-1)
	    ij[ll,:] = np.concatenate([[np.repeat(k,n_samples-k-1,0)], [np.arange(k+1,n_samples).T]]).T
	    D[ll,:] = X_[k,:] - X_[(k+1):n_samples,:]
	if (np.min(np.sum(np.abs(D),1)) == 0.) and (corr != corriid):
	    raise Exception, "Multiple X are not allowed"
	
	# Regression matrix and parameters
	F = self.regr(X)
	n_samples_F = F.shape[0]
	if F.ndim > 1:
	    p = F.shape[1]
	else:
	    p = 1
	if n_samples_F != n_samples:
	    raise Exception, "Number of rows in F and X do not match. Most likely something is wrong in the regression model."
	if p > n_samples_F:
	    raise Exception, "Ordinary least squares problem is undetermined. n_samples=%d must be greater than the regression model size p=%d!" % (n_samples, p)
	if self.beta0 != None and (self.beta0.shape[0] != p or self.beta0.ndim > 1):
	    raise Exception, "Shapes of beta0 and F do not match."
    
	# Set attributes
	self.X = X
	self.y = y
	self.theta0 = theta0
	self.thetaL = thetaL
	self.thetaU = thetaU
	self.optimizer = optimizer
	self.random_start = int(random_start)
	self.X_ = X_
	self.y_ = y_
	self.D = D
	self.ij = ij
	self.F = F
	self.nugget = nugget
	self.X_sc = np.concatenate([[mean_X],[std_X]])
	self.y_sc = np.concatenate([[mean_y],[std_y]])
	
	# Determine Kriging Model parameters
	if thetaL != None and thetaU != None:
	    # Maximum Likelihood Estimation of the parameters
	    if self.verbose:
		print "Performing Maximum Likelihood Estimation of the autocorrelation parameters..."
	    self.theta, self.score_value, self.par  = self.__arg_max_score__()
	    if np.isinf(self.score_value) :
		raise Exception, "Bad parameter region. Try increasing upper bound"
	else:
	    # Given parameters
	    if self.verbose:
		print "Given autocorrelation parameters. Computing Kriging Model parameters..."
	    self.theta = theta0
	    self.score_value, self.par = self.score()
	    if np.isinf(self.score_value):
		raise Exception, "Bad point. Try increasing theta0"
	
	if self.storage_mode == 'light':
	    # Delete heavy data (it will be computed again if required)
	    # (it is required only when MSE is wanted in self.predict)
	    if self.verbose:
		print "Light storage mode specified. Flushing autocorrelation matrix..."
	    self.D = None
	    self.ij = None
	    self.F = None
	    self.par['C'] = None
	    self.par['Ft'] = None
	    self.par['G'] = None
	
	return None
    
    
    def predict(self, X, eval_MSE = False, batch_size = None):
	"""
	This function evaluates the Kriging Model at location(s) x.
	
	Parameters
	----------
	X : array_like
	    An array with shape (n_eval, n_features) giving the point(s) at which the prediction(s) should be made.
	eval_MSE : boolean, optional
	    A boolean specifying whether the Mean Squared Error should be evaluated or not.
	    Default assumes evalMSE = False and evaluates only the BLUP (mean prediction).
	verbose : boolean, optional
	    A boolean specifying the verbose level.
	    Default is verbose = True.
	batch_size : integer, optional
	    An integer giving the maximum number of points that can be evaluated simulatneously (depending on the available memory).
	    Default is None so that all given points are evaluated at the same time.
	
	Returns
	-------
	y : array_like
	    An array with shape (n_eval, ) with the Best Linear Unbiased Prediction at x.
	MSE : array_like, optional (if eval_MSE == True)
	    An array with shape (n_eval, ) with the Mean Squared Error at x.
	"""
	
	# Check
	if np.any(np.isnan(self.par['beta'])):
	    raise Exception, "Not a valid KrigingModel. Try fitting it again with different parameters theta"
	
	# Check design & trial sites
	n_samples = self.X_.shape[0]
	if self.X_.ndim > 1:
	    n_features = self.X_.shape[1]
	else:
	    n = 1
	n_eval = X.shape[0]
	if X.ndim > 1:
	    n_features_X = X.shape[1]
	else:
	    n_features_X = 1
	X = X.reshape(n_eval, n_features_X)
	
	if n_features_X != n_features:
	    raise ValueError, "The number of features in X (X.shape[1] = %d) should match the one in self.X which is %d." % (n_features_X, n_features)
	
	if batch_size == None: # No memory management (evaluates all given points in a single batch run)
	    
	    # Normalize trial sites
	    X_ = (X - self.X_sc[0,:]) / self.X_sc[1,:]
	    
	    # Initialize output
	    y = np.zeros(n_eval)
	    if eval_MSE:
		MSE = np.zeros(n_eval)
	    
	    # Get distances to design sites
	    dx = np.zeros([n_eval*n_samples,n_features])
	    kk = np.arange(n_samples).astype(int)
	    for k in range(n_eval):
		dx[kk,:] = X_[k,:] - self.X_
		kk = kk + n_samples
	    
	    # Get regression function and correlation
	    f = self.regr(X_)
	    r = self.corr(self.theta, dx).reshape(n_eval, n_samples).T
	    
	    # Scaled predictor
	    y_ = np.matrix(f) * np.matrix(self.par['beta']) + (np.matrix(self.par['gamma']) * np.matrix(r)).T
	    
	    # Predictor
	    y = (self.y_sc[0,:] + self.y_sc[1,:] * np.array(y_)).reshape(n_eval)
	    
	    # Mean Squared Error
	    if eval_MSE:
		par = self.par
		if par['C'] == None:
		    # Light storage mode (need to recompute C, F, Ft and G)
		    if self.verbose:
			print "This KrigingModel used light storage mode at instanciation. Need to recompute autocorrelation matrix..."
		    score_value, par = self.score()
		
		rt = linalg.solve(np.matrix(par['C']), np.matrix(r))
		if self.beta0 == None:
		    # Universal Kriging
		    u = linalg.solve(np.matrix(-self.par['G'].T), np.matrix(self.par['Ft']).T * np.matrix(rt) - np.matrix(f).T)
		else:
		    # Ordinary Kriging
		    u = 0. * y
		
		MSE = self.par['sigma2'] * (1. - colsum(np.array(rt)**2.) + colsum(np.array(u)**2.)).T
		
		# Mean Squared Error might be slightly negative depending on machine precision
		# Force to zero!
		MSE[MSE < 0.] = 0.
		
		return y, MSE
		
	    else:
		
		return y
	    
	else: # Memory management
	    
	    if type(batch_size) is not int or batch_size <= 0:
		raise Exception, "batch_size must be a positive integer"
	    
	    if eval_MSE:
		
		y, MSE = np.zeros(n_eval), np.zeros(n_eval)
		for k in range(n_eval/batch_size):
		    y[k*batch_size:min([(k+1)*batch_size+1, n_eval+1])], MSE[k*batch_size:min([(k+1)*batch_size+1, n_eval+1])] = \
			self.predict(X[k*batch_size:min([(k+1)*batch_size+1, n_eval+1]),:], eval_MSE = eval_MSE, batch_size = None)
		
		return y, MSE
		
	    else:
		
		y = np.zeros(n_eval)
		for k in range(n_eval/batch_size):
		    y[k*batch_size:min([(k+1)*batch_size+1, n_eval+1])] = \
			self.__call__(x[k*batch_size:min([(k+1)*batch_size+1, n_eval+1]),:], eval_MSE = eval_MSE, batch_size = None)
		
		return y
    
    
    
    def score(self, theta = None):
	"""
	This function determines the BLUP parameters and evaluates the score function for the given autocorrelation parameters theta.
	Maximizing this score function wrt the autocorrelation parameters theta is equivalent to maximizing the likelihood of the
	assumed joint Gaussian distribution of the observations y evaluated onto the design of experiments X.
	
	Parameters
	----------
	theta : array_like, optional
	    An array containing the autocorrelation parameters at which the Kriging Model parameters should be determined.
	    Default uses the built-in autocorrelation parameters (ie theta = self.theta).
	
	Returns
	-------
	par : dict
	    A dictionary containing the requested Kriging Model parameters:
	    par['sigma2'] : Gaussian Process variance.
	    par['beta'] : Generalized least-squares regression weights for Universal Kriging or given beta0 for Ordinary Kriging.
	    par['gamma'] : Gaussian Process weights.
	    par['C'] : Cholesky decomposition of the correlation matrix [R].
	    par['Ft'] : Solution of the linear equation system : [R] x Ft = F
	    par['G'] : QR decomposition of the matrix Ft.
	    par['detR'] : Determinant of the correlation matrix raised at power 1/n_samples.
	score_value : double
	    The value of the score function associated to the given autocorrelation parameters theta.
	"""
	
	if theta == None:
	    # Use built-in autocorrelation parameters
	    theta = self.theta
	
	score_value = -np.inf
	par = { }
	
	# Retrieve data
	n_samples = self.X_.shape[0]
	D = self.D
	ij = self.ij
	F = self.F
	
	if D == None:
	    # Light storage mode (need to recompute D, ij and F)
	    if self.X_.ndim > 1:
		n_features = self.X_.shape[1]
	    else:
		n_features = 1
	    mzmax = n_samples*(n_samples-1)/2
	    ij = zeros([mzmax, n_features])
	    D = zeros([mzmax, n_features])
	    ll = array([-1])
	    for k in range(n_samples-1):
		ll = ll[-1]+1 + range(n_samples-k-1)
		ij[ll,:] = np.concatenate([[repeat(k,n_samples-k-1,0)], [np.arange(k+1,n_samples).T]]).T
		D[ll,:] = self.X_[k,:] - self.X_[(k+1):n_samples,:]
	    if (min(sum(abs(D),1)) == 0.) and (self.corr != corriid):
		raise Exception, "Multiple X are not allowed"
	    F = self.regr(self.X_)
	
	# Set up R
	r = self.corr(theta, D)
	R = np.eye(n_samples) * (1. + self.nugget)
	R[ij.astype(int)[:,0], ij.astype(int)[:,1]] = r
	R[ij.astype(int)[:,1], ij.astype(int)[:,0]] = r
	
	# Cholesky decomposition of R
	try:
	    C = linalg.cholesky(R).T
	except linalg.LinAlgError:
	    return score_value, par
	
	# Get least squares solution
	C = C.T
	Ft = linalg.solve(C, F)
	Q, G = linalg.qr(Ft)
	rcondG = 1./(linalg.norm(G)*linalg.norm(linalg.inv(G)))
	if rcondG < 1e-10:
	    # Check F
	    condF = linalg.norm(F)*linalg.norm(linalg.inv(F))
	    if condF > 1e15:
		raise Exception, "F is too ill conditioned. Poor combination of regression model and observations."
	    else:
		# Ft is too ill conditioned
		return score_value, par
	
	Yt = linalg.solve(C,self.y_)
	if self.beta0 == None:
	    # Universal Kriging
	    beta = linalg.solve(G, np.matrix(Q).T*np.matrix(Yt))
	else:
	    # Ordinary Kriging
	    beta = np.array(self.beta0)
	
	rho = np.matrix(Yt) - np.matrix(Ft) * np.matrix(beta)
	normalized_sigma2 = (np.array(rho)**2.).sum(axis=0)/n_samples
	# The determinant of R is equal to the squared product of the diagonal elements of its Cholesky decomposition C
	detR = (np.array(np.diag(C))**(2./n_samples)).prod()
	score_value = - normalized_sigma2.sum() * detR
	par = { 'sigma2':normalized_sigma2 * self.y_sc[1]**2, \
		'beta':beta, \
		'gamma':linalg.solve(C.T,rho).T, \
		'C':C, \
		'Ft':Ft, \
		'G':G }
	
	return score_value, par

    def __arg_max_score__(self):
	"""
	This function estimates the autocorrelation parameters theta as the maximizer of the score function.
	(Minimization of the opposite score function is used for convenience)
	
	Parameters
	----------
	self: All parameters are stored in the Kriging Model object.
	
	Returns
	-------
	optimal_theta : array_like
	    The best set of autocorrelation parameters (the maximizer of the score function).
	optimal_score_value : double
	    The optimal score function value.
	optimal_par : dict
	    The BLUP parameters associated to thetaOpt.
	"""
	
	if self.verbose:
	    print "The chosen optimizer is: "+str(self.optimizer)
	    if self.random_start > 1:
		print str(self.random_start)+" random starts are required."
	
	percent_completed = 0.
	
	if self.optimizer == 'fmin_cobyla':
	    
	    minus_score = lambda log10t: - self.score(theta = 10.**log10t)[0]
	    
	    constraints = []
	    for i in range(self.theta0.size):
		constraints.append(lambda log10t: log10t[i] - np.log10(self.thetaL[i]))
		constraints.append(lambda log10t: np.log10(self.thetaU[i]) - log10t[i])
	    
	    for k in range(self.random_start):
		
		if k == 0:
		    # Use specified starting point as first guess
		    theta0 = self.theta0
		else:
		    # Generate a random starting point log10-uniformly distributed between bounds
		    log10theta0 = np.log10(self.thetaL) + random.rand(self.theta0.size).reshape(self.theta0.shape) * np.log10(self.thetaU/self.thetaL)
		    theta0 = 10.**log10theta0
		
		log10_optimal_theta = optimize.fmin_cobyla(minus_score, np.log10(theta0), constraints, iprint=0)
		optimal_theta = 10.**log10_optimal_theta
		optimal_minus_score_value, optimal_par = self.score(theta = optimal_theta)
		optimal_score_value = - optimal_minus_score_value
		
		if k > 0:
		    if optimal_score_value > best_optimal_score_value:
			best_optimal_score_value = optimal_score_value
			best_optimal_par = optimal_par
			best_optimal_theta = optimal_theta
		else:
		    best_optimal_score_value = optimal_score_value
		    best_optimal_par = optimal_par
		    best_optimal_theta = optimal_theta
		if self.verbose and self.random_start > 1:
		    if (20*k)/self.random_start > percent_completed:
			percent_completed = (20*k)/self.random_start
			print str(5*percent_completed)+"% completed"
    
	else:
	    
	    raise NotImplementedError, "This optimizer ('%s') is not implemented yet. Please contribute!" % self.optimizer
    
	return best_optimal_theta, best_optimal_score_value, best_optimal_par