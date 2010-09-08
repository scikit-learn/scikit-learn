
import numpy as np
from scipy import linalg

from .base import LinearModel
from ..utils.extmath import fast_logdet


class BayesianRidge(LinearModel):
    """
    Bayesian ridge regression.

    Fit a Ridge model and optimize the regularization parameters
    lambda (precision of the weights) and alpha (precision of the
    noise) within a simple bayesian framework (MAP).

    Parameters
    ----------
    X : numpy array of shape (length,features)
        Training vectors.

    Y : numpy array of shape (length)
        Target values for training vectors.

    n_iter : int (defaut is 300)
        Stop the algorithm after a given number of steps.

    th_w : float (defaut is 1.e-12)
        Stop the algorithm if w has converged.

    compute_ll  : boolean (default is False).
        If True, compute the log-likelihood at each step of the model.

    fit_itercept : boolean, optional (True by default)

    Attributes
    ----------
    coef_ : array-like
        Mean of the weights distribution.

    lambda_ : float
        Precision parameter.

    beta_ : float
        Precision of the noise.

    sigma_ : array-like, shape [n_features,n_features]
        Variance-covariance matrix of the weights

    log_likelihood : list of float of size steps.
          Compute (if asked) the log-likelihood of the model.

    """
    def __init__(self, n_iter=300, th_w=1.e-12, compute_ll=False,
        fit_intercept=True):
        self.n_iter = n_iter
        self.th_w = th_w
        self.compute_ll = compute_ll
        self.fit_intercept = fit_intercept

    def compute_log_likelihood(self,X,Y):
        ll =  0.5 * X.shape[1] * np.log(self.lambda_)\
            + 0.5 * X.shape[0] * np.log(self.alpha_)\
            - 0.5 * self.alpha_ *  np.sum((Y - np.dot(X, self.coef_))**2)\
            - self.lambda_ * np.dot(self.coef_.T,self.coef_)\
            - 0.5 * fast_logdet(self.sigma_)\
            - 0.5 * X.shape[0] * np.log(2*np.pi)
        return ll

    def fit(self, X, Y, **params):
        """
        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data
        Y : numpy array of shape [n_samples]
            Target values

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)
        n_samples, n_features = X.shape

        X, Y = self._center_data (X, Y)

        ### "Dummy" initialization of the values of the parameters
        self.alpha_ = 1. / np.var(Y)
        self.lambda_ = 1.0
        self.log_likelihood_ = []
        U, S, V = linalg.svd(X, full_matrices=False)
        self.eigen_vals_ = S**2
        self.X_XT = np.dot(X, X.T)
        self.XT_Y = np.dot(X.T, Y)

        ### Convergence loop of the bayesian ridge regression
        for iter_ in range(self.n_iter):

            ### Compute mu and sigma (using Woodbury matrix identity)
            self.sigma_ =  np.dot(linalg.pinv(np.eye(n_samples) / self.alpha_ +
                                  self.X_XT / self.lambda_), X / self.lambda_)
            self.sigma_ = - np.dot(X.T/self.lambda_, self.sigma_)
            self.sigma_.flat[::(self.sigma_.shape[1]+1)] += 1. / self.lambda_
            self.coef_ = self.alpha_ * np.dot(self.sigma_, self.XT_Y)

            ### Update alpha and lambda
            self.gamma_ =  np.sum((self.alpha_ * self.eigen_vals_)\
                            /(self.lambda_ + self.alpha_ * self.eigen_vals_))
            self.lambda_ = self.gamma_ / np.dot(self.coef_.T,self.coef_)
            self.alpha_ = (n_samples - self.gamma_)\
                          /np.sum((Y - np.dot(X, self.coef_))**2)

            ### Compute the log likelihood
            if self.compute_ll:
                self.log_likelihood_.append(self.compute_log_likelihood(X,Y))

            ### Check for convergence
            if iter_ != 0 and np.sum(self.coef_old_ - self.coef_) < self.th_w:
                    break
            self.coef_old_ = np.copy(self.coef_)

        self.intercept_ = self._ymean - np.dot(self._xmean, self.coef_)
        # Store explained variance for __str__
        self.explained_variance_ = self._explained_variance(X, Y)
        return self




class ARDRegression(LinearModel):
    """Bayesian ard-based regression.

    Fit an ARD model and Optimize the regularization parameters alpha
    (vector of precisions of the weights) and beta (precision of the
    noise).


    Parameters
    ----------
    X : numpy array of shape (length,features)
        Training vectors.

    Y : numpy array of shape (length)
        Target values for training vectors

    n_iter : int (defaut is 300)
        Maximum number of interations.

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
    See Bishop chapter 7.2. for more details.
    This should be resived. It is not efficient and I wonder if we
    can't use libsvm for this.
    """
    # TODO: add intercept

    def __init__(self, n_iter=300, th_w=1.e-12, th_lb=1.e-12, compute_ll=False,
        fit_intercept=True):
        self.n_iter = n_iter
        self.th_w = th_w
        self.th_lb = th_lb
        self.compute_ll = compute_ll
        self.fit_intercept = fit_intercept


    def fit(self, X, Y, **params):
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)
        n_samples, n_features = X.shape

        if self.fit_intercept:
            self._xmean = X.mean(axis=0)
            self._ymean = Y.mean(axis=0)
            X = X - self._xmean
            Y = Y - self._ymean
        else:
            self._xmean = 0.
            self._ymean = 0.


        ### "Dummy" initialization of the values of the parameters
        self.alpha_ = 1./np.var(Y)
        self.lambda_ = np.ones(n_features)
        self.log_likelihood_ = []
        self.X_XT = np.dot(X,X.T)
        self.XT_Y = np.dot(X.T,Y)

        ### Launch the convergence loop
        self.loop_ard(X,Y)
        
        self.intercept_ = self._ymean - np.dot(self._xmean, self.coef_)
        # Store explained variance for __str__
        self.explained_variance_ = self._explained_variance(X, Y)
        return self

    
    def loop_ard(self,X,Y):
      
        n_samples, n_features = X.shape
      
        ### Convergence loop of the bayesian ridge regression
        for iter_ in range(self.n_iter):

            ### Compute mu and sigma (using Woodbury matrix identity)
            self.sigma_ =  np.dot(
                    linalg.pinv(np.eye(n_samples)/self.alpha_ +
                    np.dot(X * np.reshape(1./self.lambda_,[1,-1]),X.T)),
                    X * np.reshape(1./self.lambda_,[1,-1]))
            self.sigma_ = - np.dot(np.reshape(1./self.lambda_,[-1,1]) * X.T ,
                                                                    self.sigma_)
            self.sigma_.flat[::(self.sigma_.shape[1]+1)] += 1./self.lambda_
            self.coef_ = self.alpha_ *np.dot(self.sigma_,self.XT_Y)

            ### Update alpha and lambda
            self.gamma_ =  1. - self.lambda_*np.diag(self.sigma_)
            coef_2 = (self.coef_)**2
            self.lambda_[coef_2>self.th_lb] = self.gamma_[coef_2>self.th_lb]\
                                                   /coef_2[coef_2>self.th_lb]
            self.lambda_[coef_2 <= self.th_lb] = 1./self.th_lb
            self.alpha_ = (n_samples - self.gamma_.sum())\
                          /np.sum((Y - np.dot(X, self.coef_))**2)

            #### Compute the log likelihood
            if self.compute_ll:
                self.log_likelihood_.append(self.compute_log_likelihood(X,Y))

            ### Check for convergence
            if iter_ != 0 and np.sum(self.coef_old_ - self.coef_) < self.th_w:
                    break
            self.coef_old_ = np.copy(self.coef_)




 #### Compute the log likelihood
            #if ll_bool :
                #A_ = np.eye(X.shape[1])/alpha
                #C_ = (1./beta)*np.eye(X.shape[0]) + np.dot(X,np.dot(A_,X.T))
                #ll = X.shape[0]*np.log(2*np.pi)+fast_logdet(C_)
                #ll += np.dot(Y.T,np.dot(linalg.pinv(C_),Y))
                #log_likelihood.append(-0.5*ll)


class RVM_R(ARDRegression):
    """
    See Bishop chapter 7.2. for more details.
    This should be resived. It is not efficient and I wonder if we
    can't use libsvm for this.
    """

    def __init__(self, n_iter=300, th_w=1.e-12, th_lb=1.e-12, kernel = "linear",
                        compute_ll=False, fit_intercept=True):
        self.n_iter = n_iter
        self.th_w = th_w
        self.th_lb = th_lb
        self.kernel = kernel
        self.compute_ll = compute_ll
        self.fit_intercept = fit_intercept


    def fit(self, X, Y, **params):
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)
        if self.kernel == "linear":
            self.X_save = X
            X = np.dot(X,X.T)
        n_samples, n_features = X.shape

        if self.fit_intercept:
            self._xmean = X.mean(axis=0)
            self._ymean = Y.mean(axis=0)
            X = X - self._xmean
            Y = Y - self._ymean
        else:
            self._xmean = 0.
            self._ymean = 0.


        ### "Dummy" initialization of the values of the parameters
        self.alpha_ = 1./np.var(Y)
        self.lambda_ = np.ones(n_features)
        self.log_likelihood_ = []
        self.X_XT = np.dot(X,X.T)
        self.XT_Y = np.dot(X.T,Y)

        ### Launch the convergence loop
        self.loop_ard(X,Y)
        
        self.intercept_ = self._ymean - np.dot(self._xmean, self.coef_)
        # Store explained variance for __str__
        #self.explained_variance_ = self._explained_variance(X, Y)
        return self

    def predict(self,Xtest):
        Xtest = np.asanyarray(Xtest)
        Xtest = np.dot(self.X_save,Xtest.T)
        return np.dot(Xtest, self.coef_) + self.intercept_
        


