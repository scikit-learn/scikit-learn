
import numpy as np
from scipy import optimize
from abc import abstractmethod
from sklearn.linear_model.base import LinearModel
from .shapley_fast import ShapleyValue, wcov, wcorr

class ConstrainedRegression(LinearModel):
    """
    Base class for constrained regression models.
    """

    def __init__(self, fit_intercept = True):
        self.fit_intercept = fit_intercept
        self.coef_ = None
        self.intercept_ = 0
        self.importances = None

    def constrained_optimization(self, corr):
        """
        Find the linear regression coefficients given desired net effects. A wrapper for a bounded L-BFGS-B
        optimizer.
        """
        fit = optimize.minimize( lambda x: np.sum( (np.multiply(x.dot(corr[1:,1:]), x) - self.importances)**2 ), 
                    method='L-BFGS-B',
                    x0 = np.array([0.001]*corr[1:,1:].shape[0]), 
                    bounds = [(0,None)]*corr[1:,1:].shape[0] )
        return fit

    @abstractmethod
    def fit(self, X, y, weights = None ):
        """
        Fit  model
        """

    def decision_function(self, X):
        """
        Decision function of the linear model
        """
        return self.predict(X)

    def predict(self, X):
        """
        Predict using the linear model
        """
        return self.intercept_ + X.dot(self.coef_)


    def score(self,X, y, weights = None):
        """
        Returns the coefficient of determination R^2 of the prediction.
        """

        if weights is None: weights = np.ones(y.shape[0])

        y_mean = weights/np.sum(weights) *y
        y_pred = self.predict(X)
        u = ((y - y_pred)**2).sum()
        v = ((y - y_mean)**2).sum()
        return 1.-u/v


class ShapleyRegression(ConstrainedRegression):
    """
    Shapley valued regression
    """

    def fit(self, X, y, weights = None ):
        """
        Fit a linear regression with prescribed importances defined by the Shapley value for each covariate
        """
        if weights is None: weights = np.ones(y.shape[0])
        data = np.hstack((y.reshape(y.shape[0],1),X))

        S = wcov(data, weights)
        corr = wcorr(data, weights)
        wsd = np.sqrt(S.diagonal())

        self.importances = ShapleyValue( S )
        model = self.constrained_optimization( corr )

        if self.fit_intercept:
            w = np.diagflat( weights/np.sum(weights),k=0)
            wmean = np.sum(w.dot(data), axis=0)
            self.intercept_ = wmean[0] - wsd[0]*np.sum(wmean[1:]*model.x/wsd[1:])

        self.coef_ = wsd[0]*model.x/wsd[1:] 

        return self