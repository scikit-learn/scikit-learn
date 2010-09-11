"""
Various bayesian regression
"""

# Authors: V. Michel, F. Pedregosa
# License: BSD 3 clause


import numpy as np
from scipy import linalg

from .base import LinearModel
from ..utils.extmath import fast_logdet

###############################################################################
# BayesianRidge regression

class BayesianRidge(LinearModel):
    """Bayesian ridge regression

    Fit a Bayesian ridge model and optimize the regularization parameters
    lambda (precision of the weights) and alpha (precision of the noise).

    Parameters
    ----------
    X : numpy array of shape (length,features)
        Training vectors.

    Y : numpy array of shape (length)
        Target values for training vectors

    n_iter : int (default is 300)
        Maximum number of interations.

    eps : float (default is 1.e-3)
        Stop the algorithm if w has converged.

    alpha_1 : float (default is 1.e-6)
        Hyper-parameter : shape parameter for the Gamma distribution prior over
        the alpha parameter.

    alpha_2 : float (default is 1.e-6)
        Hyper-parameter : inverse scale parameter (rate parameter) for the Gamma
        distribution prior over the alpha parameter.

    lambda_1 : float (default is 1.e-6)
        Hyper-parameter : shape parameter for the Gamma distribution prior over
        the lambda parameter.

    lambda_2 : float (default is 1.e-6)
        Hyper-parameter : inverse scale parameter (rate parameter) for the Gamma
        distribution prior over the lambda parameter.

    compute_score : boolean (default is False)
        If True, compute the objective function at each step of the model.

    fit_intercept : boolean (default is True)
        wether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    Attributes
    ----------
    coef_ : numpy array of shape (nb_features)
        Coefficients of the regression model (mean of the weights
        distribution.)

    alpha_ : float
       estimated precision of the noise.

    lambda_ : numpy array of shape (nb_features)
       estimated precisions of the weights.

    score_ : float
        if computed, value of the objective function (to be maximized)

    Methods
    -------
    fit(X, y) : self
        Fit the model

    predict(X) : array
        Predict using the model.


    Examples
    --------

    """


    def __init__(self, n_iter=300, eps=1.e-3, alpha_1 = 1.e-6, alpha_2 = 1.e-6,
                lambda_1=1.e-6, lambda_2=1.e-6, compute_score=False,
                fit_intercept=True):
        """
        Parameters
        ----------
        n_iter : int (default is 300)
            Maximum number of interations.

        eps : float (default is 1.e-3)
            Stop the algorithm if w has converged.

        alpha_1 : float (default is 1.e-6)
            Hyper-parameter : shape parameter for the Gamma distribution prior
            over the alpha parameter.

        alpha_2 : float (default is 1.e-6)
            Hyper-parameter : inverse scale parameter (rate parameter) for the
            Gamma distribution prior over the alpha parameter.

        lambda_1 : float (default is 1.e-6)
            Hyper-parameter : shape parameter for the Gamma distribution prior
            over the lambda parameter.

        lambda_2 : float (default is 1.e-6)
            Hyper-parameter : inverse scale parameter (rate parameter) for the
            Gamma distribution prior over the lambda parameter.

        compute_score : boolean (default is False)
            If True, compute the objective function at each step of the model.

        fit_intercept : boolean (default is True)
          wether to calculate the intercept for this model. If set
          to false, no intercept will be used in calculations
          (e.g. data is expected to be already centered).

        """
        self.n_iter = n_iter
        self.eps = eps
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.compute_score = compute_score
        self.fit_intercept = fit_intercept

    def fit(self, X, Y, verbose=False, **params):
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
        X, Y, Xmean, Ymean = self._center_data (X, Y)
        n_samples, n_features = X.shape

        ### Initialization of the values of the parameters
        self.alpha_ = 1./np.var(Y)
        self.lambda_ = 1.
        XT_Y = np.dot(X.T, Y)
        U, S, Vh = linalg.svd(X, full_matrices=False)
        eigen_vals_ = S**2
        self.all_score_ = []

        ### Convergence loop of the bayesian ridge regression
        for iter_ in range(self.n_iter):

            ### Compute mu and sigma
            # sigma_ = lambda_ / alpha_ * np.eye(n_features) + np.dot(X.T, X)
            # coef_ = sigma_^-1 * XT * y
            if n_samples > n_features:
                coef_ = np.dot(Vh.T,
                    Vh / (eigen_vals_ + self.lambda_ / self.alpha_)[:,None])
                coef_ = np.dot(coef_, XT_Y)
                if self.compute_score:
                    logdet_sigma_ = - np.sum(
                        np.log(self.lambda_ + self.alpha_* eigen_vals_))
            else:
                coef_ = np.dot(X.T, np.dot(U /
                    (eigen_vals_ + self.lambda_ / self.alpha_)[None,:],
                    U.T))
                coef_ = np.dot(coef_, Y)
                if self.compute_score:
                    logdet_sigma_ = self.lambda_ * np.ones(n_features)
                    logdet_sigma_[:n_samples] += self.alpha_ * eigen_vals_
                    logdet_sigma_ = - np.sum(np.log(logdet_sigma_))

            if self.compute_score:
                self.logdet_sigma_ = logdet_sigma_

            ### Update alpha and lambda
            self.rmse_ = np.sum((Y - np.dot(X, coef_))**2)
            self.gamma_ =  np.sum((self.alpha_ * eigen_vals_) \
                            / (self.lambda_ + self.alpha_ * eigen_vals_))
            self.lambda_ =  (self.gamma_ + 2*self.lambda_1) \
                            / (np.sum(coef_**2) + 2*self.lambda_2)
            self.alpha_ = (n_samples - self.gamma_ +  2*self.alpha_1) \
                          / (self.rmse_ + 2*self.alpha_2)

            self.coef_ = coef_

            ### Compute the objective function
            if self.compute_score:
                self.all_score_.append(self.objective_function(X))

            ### Check for convergence
            if (iter_ != 0 and np.sum(np.abs(coef_old_ - coef_)) < self.eps):
                if verbose:
                    print "Convergence after ", str(iter_), " iterations"
                break
            coef_old_ = np.copy(coef_)

        self._set_intercept(Xmean, Ymean)
        # Store explained variance for __str__
        self.explained_variance_ = self._explained_variance(X, Y)
        return self


    def objective_function(self, X):
        """
        Compute the objective function.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        Returns
        -------
        score_ : value of the objective function (to be maximized)
        """
        score_ = self.lambda_1 * np.log(self.lambda_) - self.lambda_2 \
                       * self.lambda_
        score_ += self.alpha_1 * np.log(self.alpha_) - self.alpha_2 \
                       * self.alpha_
        score_ += 0.5 * X.shape[1] * np.log(self.lambda_) \
                       + 0.5 * X.shape[0] * np.log(self.alpha_) \
                       - 0.5 * self.alpha_ *  self.rmse_ \
                       - 0.5 * (self.lambda_ * np.sum(self.coef_**2)) \
                       - 0.5 * self.logdet_sigma_ \
                       - 0.5 * X.shape[0] * np.log(2*np.pi)
        return score_


###############################################################################
# ARD (Automatic Relevance Determination) regression


class ARDRegression(LinearModel):
    """Bayesian ARD regression.

    Fit the weights of a regression model, using an ARD prior. The weights of
    the regression model are assumed to be in Gaussian distributions.
    Also estimate the parameters lambda (precisions of the distributions of the
    weights) and alpha (precision of the distribution of the noise).
    The estimation is done by an iterative procedures (Evidence Maximization)

    Parameters
    ----------
    X : numpy array of shape (length,features)
        Training vectors.

    Y : numpy array of shape (length)
        Target values for training vectors

    n_iter : int (default is 300)
        Maximum number of interations.

    eps : float (default is 1.e-3)
        Stop the algorithm if w has converged.

    alpha_1 : float (default is 1.e-6)
        Hyper-parameter : shape parameter for the Gamma distribution prior over
        the alpha parameter.

    alpha_2 : float (default is 1.e-6)
        Hyper-parameter : inverse scale parameter (rate parameter) for the Gamma
        distribution prior over the alpha parameter.

    lambda_1 : float (default is 1.e-6)
        Hyper-parameter : shape parameter for the Gamma distribution prior over
        the lambda parameter.

    lambda_2 : float (default is 1.e-6)
        Hyper-parameter : inverse scale parameter (rate parameter) for the Gamma
        distribution prior over the lambda parameter.

    compute_score : boolean (default is False)
        If True, compute the objective function at each step of the model.

    threshold_lambda : float (default is 1.e+4)
        threshold for removing (pruning) weights with high precision from
        the computation.

    fit_intercept : boolean (default is True)
        wether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    Attributes
    ----------
    coef_ : numpy array of shape (nb_features)
        Coefficients of the regression model (mean of the weights
        distribution.)

    alpha_ : float
       estimated precision of the noise.

    lambda_ : numpy array of shape (nb_features)
       estimated precisions of the weights.

    sigma_ : numpy array of shape (nb_features,nb_features)
        estimated variance-covariance matrix of the weights

    score_ : float
        if computed, value of the objective function (to be maximized)

    Methods
    -------
    fit(X, y) : self
        Fit the model

    predict(X) : array
        Predict using the model.


    Examples
    --------
    """

    def __init__(self, n_iter=300, eps=1.e-3, alpha_1 = 1.e-6, alpha_2 = 1.e-6,
                  lambda_1 = 1.e-6, lambda_2 = 1.e-6, compute_score = False,
                  threshold_lambda = 1.e+4, fit_intercept = True):
        """
        Parameters
        ----------
        n_iter : int (default is 300)
            Maximum number of interations.

        eps : float (default is 1.e-3)
            Stop the algorithm if w has converged.

        alpha_1 : float (default is 1.e-6)
            Hyper-parameter : shape parameter for the Gamma distribution prior
            over the alpha parameter.

        alpha_2 : float (default is 1.e-6)
            Hyper-parameter : inverse scale parameter (rate parameter) for the
            Gamma distribution prior over the alpha parameter.

        lambda_1 : float (default is 1.e-6)
            Hyper-parameter : shape parameter for the Gamma distribution prior
            over the lambda parameter.

        lambda_2 : float (default is 1.e-6)
            Hyper-parameter : inverse scale parameter (rate parameter) for the
            Gamma distribution prior over the lambda parameter.

        compute_score : boolean (default is False)
            If True, compute the objective function at each step of the model.

        threshold_lambda : float (default is 1.e+4)
            threshold for removing (pruning) weights with high precision from
            the computation.

        fit_intercept : boolean (default is True)
            wether to calculate the intercept for this model. If set
            to false, no intercept will be used in calculations
            (e.g. data is expected to be already centered).
        """
        self.n_iter = n_iter
        self.eps = eps
        self.fit_intercept = fit_intercept
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.compute_score = compute_score
        self.threshold_lambda = threshold_lambda


    def fit(self, X, Y, **params):
        """
        Fit the ARDRegression model according to the given training data and
        parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        Y : array, shape = [n_samples]
            Target values (integers)

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)

        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)

        n_samples, n_features = X.shape

        X, Y, Xmean, Ymean = self._center_data (X, Y)

        ### Initialization of the values of the parameters
        self.alpha_ = 1./np.var(Y)
        self.lambda_ = np.ones(n_features)
        self.all_score_ = []

        ### Launch the convergence loop
        self.evidence_maximization(X, Y)

        self._set_intercept(Xmean, Ymean)
        # Store explained variance for __str__
        self.explained_variance_ = self._explained_variance(X, Y)
        return self

    def evidence_maximization(self, X, Y, verbose=False):
        """
        Iterative procedure for estimating the ARDRegression model according to
        the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        Y : array, shape = [n_samples]
            Target values (integers)

        Attributes
        ----------
        keep_lambda : boolean numpy array of shape (nb_features)
            Lambda under a given threshold, to be keep for the computation.
            Avoid divergence when lambda is to high
        """

        n_samples, n_features = X.shape
        coef_ = np.zeros(n_features)
        keep_lambda = np.ones(n_features,dtype=bool)

        ### Iterative procedure of ARDRegression
        for iter_ in range(self.n_iter):

            ### Compute mu and sigma (using Woodbury matrix identity)
            self.sigma_ = linalg.pinv(np.eye(n_samples)/self.alpha_ +
                          np.dot(X[:,keep_lambda] *
                          np.reshape(1./self.lambda_[keep_lambda],[1,-1]),
                          X[:,keep_lambda].T))
            self.sigma_ = np.dot(self.sigma_,X[:,keep_lambda]
                          * np.reshape(1./self.lambda_[keep_lambda],
                          [1,-1]))
            self.sigma_ = - np.dot(np.reshape(1./self.lambda_[keep_lambda],
                          [-1,1]) * X[:,keep_lambda].T ,self.sigma_)
            self.sigma_.flat[::(self.sigma_.shape[1]+1)] += \
                          1./self.lambda_[keep_lambda]
            coef_[keep_lambda] = self.alpha_ \
                            * np.dot(self.sigma_,np.dot(X[:,keep_lambda].T,
                            Y))

            ### Update alpha and lambda
            self.rmse_ = np.sum((Y - np.dot(X, coef_))**2)
            self.gamma_ =  1. - self.lambda_[keep_lambda]\
                                          *np.diag(self.sigma_)
            self.lambda_[keep_lambda] = (self.gamma_ + 2*self.lambda_1)\
                        /((coef_[keep_lambda])**2 + 2*self.lambda_2)
            self.alpha_ = (n_samples - self.gamma_.sum() +  2*self.alpha_1)\
                            /(self.rmse_ + 2*self.alpha_2)

            ### Prune the weights with a precision over a threshold
            keep_lambda = self.lambda_ < self.threshold_lambda
            coef_[keep_lambda == False] = 0

            self.coef_ = coef_

            ### Compute the objective function
            if self.compute_score:
                self.all_score_.append(self.objective_function(X))

            ### Check for convergence
            if iter_ != 0 and np.sum(np.abs(coef_old_ - coef_)) < self.eps:
                if verbose:
                    print "Convergence after %s iterations" % iter_
                break
            coef_old_ = np.copy(coef_)


    def objective_function(self, X):
        """
        Compute the objective function.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        Returns
        -------
        score_ : value of the objective function (to be maximized)
        """

        score_ = (self.lambda_1 * np.log(self.lambda_) - self.lambda_2\
                       * self.lambda_).sum()
        score_ += self.alpha_1 * np.log(self.alpha_) - self.alpha_2\
                       * self.alpha_
        score_ += 0.5 * (fast_logdet(self.sigma_)  + X.shape[0]\
                          * np.log(self.alpha_) + np.sum(np.log(self.lambda_)))
        score_ -= 0.5 * (self.alpha_ * self.rmse_\
                              + (self.lambda_ * self.coef_**2).sum())
        return score_



