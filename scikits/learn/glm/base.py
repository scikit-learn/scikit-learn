# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Vincent Michel <vincent.michel@inria.fr>
#
# License: BSD Style.


"""
Generalized Linear models.
"""

import warnings

from math import log10
import numpy as np
import scipy.linalg 

from ..utils.extmath import fast_logdet
from ..base import BaseEstimator, RegressorMixin

###
### TODO: intercept for all models
### We should define a common function to center data instead of
### repeating the same code inside each fit method.
###
### Also, bayesian_ridge_regression and bayesian_regression_ard
### should be squashed into its respective objects.
###

class LinearModel(BaseEstimator, RegressorMixin):
    """Base class for Linear Models"""

    def predict(self, X):
        """
        Predict using the linear model

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Returns predicted values.
        """
        X = np.asanyarray(X)
        return np.dot(X, self.coef_) + self.intercept_

    def _explained_variance(self, X, Y):
        """Compute explained variance a.k.a. r^2"""
        ## TODO: this should have a tests.
        return 1 - np.linalg.norm(Y - self.predict(X))**2 \
                         / np.linalg.norm(Y)**2

    def __str__(self):
        if self.coef_ is not None:
            return ("%s \n%s #... Fitted: explained variance=%s" %
                    (repr(self), ' '*len(self.__class__.__name__),  
                     self.explained_variance_))
        else:
            return "%s \n#... Not fitted to data" % repr(self)


class LinearRegression(LinearModel):
    """
    Ordinary least squares Linear Regression.

    Attributes
    ----------
    `coef_` : array
        Estimated coefficients for the linear regression problem.

    `intercept_` : array
        Independent term in the linear model.

    Notes
    -----
    From the implementation point of view, this is just plain Ordinary
    Least Squares (numpy.linalg.lstsq) wrapped as a predictor object.

    """

    def __init__(self, fit_intercept=True):
        self.fit_intercept = fit_intercept


    def fit(self, X, Y, **params):
        """
        Fit linear model.

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data
        Y : numpy array of shape [n_samples]
            Target values
        fit_intercept : boolean, optional
            wether to calculate the intercept for this model. If set
            to false, no intercept will be used in calculations
            (e.g. data is expected to be already centered).

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)
        X = np.asanyarray( X )
        Y = np.asanyarray( Y )

        if self.fit_intercept:
            # augmented X array to store the intercept
            X = np.c_[X, np.ones(X.shape[0])]
        self.coef_, self.residues_, self.rank_, self.singular_ = \
                np.linalg.lstsq(X, Y)
        if self.fit_intercept:
            self.intercept_ = self.coef_[-1]
            self.coef_ = self.coef_[:-1]
        else:
            self.intercept_ = np.zeros(self.coef_X.shape[1])
        return self

