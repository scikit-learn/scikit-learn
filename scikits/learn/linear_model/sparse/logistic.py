"""
Sparse Logistic Regression module

This module has the same API as scikits.learn.linear_model.logistic, but is
designed to handle efficiently data in sparse matrix format.
"""

import numpy as np

from ...base import ClassifierMixin
from ...svm.sparse.base import SparseBaseLibLinear
from ...linear_model.sparse.base import CoefSelectTransformerMixin
from ...svm.liblinear import csr_predict_prob

class LogisticRegression(SparseBaseLibLinear, ClassifierMixin,
                         CoefSelectTransformerMixin):
    """
    Logistic Regression.

    Implements L1 and L2 regularized logistic regression.

    Parameters
    ----------

    penalty : string, 'l1' or 'l2'
        Used to specify the norm used in the penalization

    dual : boolean
        Dual or primal formulation. Dual formulation is only
        implemented for l2 penalty.

    C : float
        Specifies the strength of the regularization. The smaller it is
        the bigger in the regularization.

    fit_intercept : bool, default: True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added the decision function

    intercept_scaling : float, default: 1
        when self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased

    tol: float, optional
         tolerance for stopping criteria

    Attributes
    ----------

    `coef_` : array, shape = [n_classes-1, n_features]
        Coefficient of the features in the decision function.

    `intercept_` : array, shape = [n_classes-1]
        intercept (a.k.a. bias) added to the decision function.
        It is available only when parameter intercept is set to True

    See also
    --------
    LinearSVC

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon,
    to have slightly different results for the same input data. If
    that happens, try with a smaller tol parameter.

    References
    ----------
    LIBLINEAR -- A Library for Large Linear Classification
    http://www.csie.ntu.edu.tw/~cjlin/liblinear/
    """

    def __init__(self, penalty='l2', dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1):

        super(LogisticRegression, self).__init__ (penalty=penalty,
            dual=dual, loss='lr', tol=tol, C=C,
            fit_intercept=fit_intercept, intercept_scaling=intercept_scaling)

    def predict_proba(self, X):
        """
        Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.
        """
        import scipy.sparse
        X = scipy.sparse.csr_matrix(X)
        X.data = np.asanyarray(X.data, dtype=np.float64, order='C')
        probas = csr_predict_prob(X.shape[1], X.data, X.indices,
                                  X.indptr, self.raw_coef_,
                                  self._get_solver_type(),
                                  self.tol, self.C,
                                  self.class_weight_label,
                                  self.class_weight, self.label_,
                                  self._get_bias())
        return probas[:,np.argsort(self.label_)]

    def predict_log_proba(self, T):
        """
        Log of Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.
        """
        return np.log(self.predict_proba(T))

    def min_C(self, X, y):
        """
        Return the maximum value for C that yields a model with coefficients
        and intercept set to zero.

        This value is valid if class_weight parameter in fit() is not set.
        min_C() raises ValueError if penalty is different from 'l1' or
        the number of classes in y is different from 2.

        Parameters
        ----------
        X : sparse matrix, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array, shape = [n_samples]
            Target vector relative to X

        Returns
        -------
        min_C: float
            minimum value for C
        """
        from scipy import sparse

        if self.penalty != 'l1':
            raise ValueError('penalty is not l1')
        classes = np.unique(y)
        if len(classes) != 2:
            raise ValueError('min_C: number of classes != 2')

        X = sparse.csc_matrix(X)
        y = np.asanyarray(y)

        if self.fit_intercept:
            bias = self.intercept_scaling * np.ones((np.size(y), 1))
            X = sparse.hstack((X, bias))

        _y = np.empty(y.shape)
        _y[y == classes[0]] = -1
        _y[y == classes[1]] = 1
        _y = _y.reshape((1, -1))

        den = np.max(np.abs(_y * X))
        if den == 0.0:
            raise ValueError('Ill-posed min_C calculation')
        return 2.0 / den
