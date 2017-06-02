"""
Algorithms for Boosting:
- Functional Gradient Descent
"""

# Authors: James Bergstra
# License: BSD3

import numpy as np

from .utils import safe_asarray
from .ensemble import BaseEnsemble


class FitIter(object):
    """
    Iterations (self.next()) implement one round of functional gradient
    boosting.

    Attributes
    ----------
    fgb : the FunctionalGradientBoosting instance
        FitIter implements the self.fit of this object.

    X : array-like of shape = [n_samples, n_features]
        Training input samples

    residual : array of shape = [n_samples]
        Running regression target (originally the training target)
    
    N.B. This object works in-place on self.residual

    """
    def __init__(self, fgb, X, residual):
        self.fgb = fgb
        self.X = X
        self.residual = residual

    def __iter__(self):
        return self

    def next(self):
        if self.fgb.n_estimators == len(self.fgb.estimators_):
            raise StopIteration
        if self.fgb.estimators_:
            self.residual -= self.fgb.estimators_[-1].predict(self.X)
        base = self.fgb._make_estimator()
        base.fit(self.X, self.residual)
        return self


class FunctionalGradientBoosting(BaseEnsemble):
    """
    Regression Boosting via functional gradient descent.

    The algorithm is to construct a regression ensemble by using a "base
    estimator" to repeatedly fit residual training error. So for example, the
    first iteration fits some function f() to the original (X, y) training
    data, and the second iteration fits some g() to (X, y - f(X)), the third
    iterations fits some h() to (X y - f(X) - g(X)), and so on.  The final
    ensemble is f() + g() + h() + ...

    This procedure is equivalent to functional gradient descent when the the
    training objective is to minimize mean squared error (MSE).

    For more information see e.g.:
    J. H. Friedman (2002). "Stochastic Gradient Boosting",
    Computational Statistics & Data Analysis.

    TODO: Mason has a good paper on the subject as well.
    """

    def __init__(self, base_estimator, n_estimators,):
        super(FunctionalGradientBoosting, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators)

    def fit_iter(self, X, y):
        """Create a fitting iterator for training set X, y.

        See class FitIter().
        """
        X = safe_asarray(X)
        y = np.array(y)    # N.B. makes a copy
        if 'int' in str(y.dtype):
            raise NotImplementedError('ints typically mean classif')
        return FitIter(self, X, y)

    def fit(self, X, y):
        """Build a regression ensemble by funtional gradient boosting.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        Return
        ------
        self : object
            Returns self.
        """
        for _ in self.fit_iter(X, y):
            pass
        return self

    def predict(self, X):
        """Return the prediction for array-like X.
        
        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            Test samples.

        Return
        ------
        prediction : numpy array of shape = [n_samples]
            Test predictions.

        """
        rval = self.estimators_[0].predict(X)
        for estimator in self.estimators_[1:]:
            pred_i = estimator.predict(X) 
            rval += pred_i
        return rval
