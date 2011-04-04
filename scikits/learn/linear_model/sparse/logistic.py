"""
Sparse Logistic Regression module

This module has the same API as scikits.learn.glm.logistic, but is
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
