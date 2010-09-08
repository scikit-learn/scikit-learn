
import numpy as np
from scipy import linalg

from .base import LinearModel
from ..utils.extmath import fast_logdet


class BayesianRidge(LinearModel):
    """
    Bayesian ridge regression. Optimize the regularization parameters lambda
    (precision of the weights) and alpha (precision of the noise) within a
    simple bayesian framework (MAP).

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

    Notes
    -----
    See Bishop p 167-169 for more details.
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
        



class VariationalBayes(object):

    def __init__(self,Q_,nb_iter_,save=True,verbose=True):
        self.Q_ = Q_
        self.nb_iter_ = nb_iter_
        self.save = save
        self.verbose = verbose
        self.priors = False

    def set_priors(self,alpha1_,alpha2_,lambda1_,lambda2_,
                        pjq_,etaq_):
        self.alpha1_ = alpha1_
        self.alpha2_ = alpha2_
        self.lambda1_ = lambda1_
        self.lambda2_ = lambda2_
        self.pjq_ = pjq_
        self.etaq_ = etaq_
        self.all_time_ = []
        self.all_coef_ = []
        self.all_pjq_ = []
        self.all_qzjq_ = []
        self.all_alpha_ = []
        self.all_lambda_ = []
        self.all_training_mse_ = []
        self.all_classes_ = []

class Variational_MVM_R(MVM):

    def predict(self,y_t):
        y_t = y_t - self.mean_X_
        self.prediction = np.dot(y_t, self.coef_) + self.ymean_
        return np.ravel(self.prediction)

    def fit(self, X, y, qzjq_=None, fixed_z = -1):

        ### Add the interecpt
        self.ymean_ = np.mean(y)
        y = y - self.ymean_
        self.y = np.reshape(y,[np.size(y),1])
        self.mean_X_ = np.mean(X,0)
        self.X = X - self.mean_X_
        self.XT_ = np.transpose(self.X)
        self.Gram_ = np.dot(self.XT_,self.X)
        self.XTy_ = np.dot(self.XT_,self.y)
        self.n_ = np.size(self.y)
        self.p_ = np.size(self.X,1)


        ### Initiate the parameters
        print "Initialize parameters"
        self.a1_ = np.copy(self.alpha1_)
        self.a2_ = np.copy(self.alpha2_)
        if qzjq_ == None :
            self.qzjq_ = (1./self.Q_)*np.ones([self.p_,self.Q_])
        else :
            self.qzjq_ = qzjq_
        self.l1_ = np.copy(self.lambda1_)
        self.l2_ = np.copy(self.lambda2_)
        self.dq_ = np.copy(self.etaq_)

        ######
        ### Estimation Loop
        ######
        for ite in range(self.nb_iter_):

            start_time_ = time.time()

            ### w
            self.A = np.sum((self.l1_/self.l2_)*self.qzjq_,1)
            invA = np.diag(1./self.A)
            invMat = nl.pinv(np.eye(self.n_)*self.a2_/self.a1_\
                      + np.dot(self.X, np.dot(invA, self.XT_)))
            invAX = np.dot(invA, self.XT_)
            self.Sigma_ = invA - np.dot(invAX, np.dot(invMat,invAX.T))
            self.mu_ = (self.a1_/self.a2_*np.dot(self.Sigma_,self.XTy_))
            self.muj_sigmajj_ = np.ravel(self.mu_)**2+np.diag(self.Sigma_)
            self.muj_sigmajj_ = np.reshape(self.muj_sigmajj_, [self.p_,1])
            self.coef_ = np.copy(self.mu_)

            ### lambda
            self.l1_ = self.lambda1_ + 0.5*np.sum(self.qzjq_, 0)
            self.l2_ = self.lambda2_ + 0.5*np.sum(self.muj_sigmajj_ *
                                                            self.qzjq_, 0)


            ### alpha
            self.a1_ = self.alpha1_ + 0.5*self.n_
            yXw = self.y - np.dot(self.X,self.mu_)
            trace = np.sum(np.multiply(self.Sigma_,self.Gram_))
            self.a2_ = np.float(self.alpha2_ + 0.5*np.dot(yXw.T, yXw) + trace)
            lambda_time_ = time.time()


            ### q(zj = q)
            if ite >= fixed_z:
                self.qzjq_ = np.exp(-0.5*self.muj_sigmajj_*self.l1_/self.l2_\
                              + np.log(self.pjq_)\
                              + 0.5*(sps.digamma(self.l1_)-np.log(self.l2_)))
                norm = np.reshape(np.sum(self.qzjq_,1),[self.p_])
                for q in range(self.Q_) :
                    self.qzjq_[norm!=0,q] /= norm[norm!=0]
                    self.qzjq_[norm==0,q] = 1.*np.ones(np.sum(norm==0))/self.Q_

                ### deltaq
                self.dq_ = self.etaq_ + np.sum(self.qzjq_,0)
                self.pjq_ = np.exp(sps.digamma(self.dq_)-
                            sps.digamma(np.sum(self.dq_)))



            #### output and save
            mse_ = np.sqrt(np.sum(yXw**2))
            self.all_time_.append(time.time() - start_time_)
            classes_ = np.copy(np.zeros(self.Q_))
            #for q in range(self.Q_):
                #classes_[q] = \
                #np.size(np.where(np.argmax(self.qzjq_,1)==q))
            #self.all_classes_.append(classes_)
            print 50*"#"
            print "Iteration : ",ite
            print "Total step time : ",self.all_time_[-1]
            print "Training MSE : ",np.sqrt(mse_)
            self.all_training_mse_.append(mse_)
            if self.save == True :
                self.all_coef_.append(np.copy(self.coef_))
                self.all_alpha_.append(self.a1_/self.a2_)
                self.all_lambda_.append(self.l1_/self.l2_)
                self.all_pjq_.append(np.copy(self.pjq_))
                self.all_qzjq_.append(np.copy(self.qzjq_))
            if self.verbose == True :
                #print "Classes : ",classes_
                print "alpha : ",self.a1_/self.a2_
                print "a1 : ",self.a1_
                print "a2 : ",self.a2_
                print "lambda : ",self.l1_/self.l2_
                print "l1 : ",self.l1_
                print "l2 : ",self.l2_
                print "dq : ",self.dq_
        return self












#def fit(self, X):
        #"""
        #Detects the soft boundary (aka soft boundary) of the set of samples X.

        #Parameters
        #----------
        #X : array-like, shape = [n_samples, n_features]
            #Set of samples, where n_samples is the number of samples and
            #n_features is the number of features.

        #"""
        #super(OneClassSVM, self).fit(X, [])

