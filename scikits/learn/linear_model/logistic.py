import numpy as np

from ..base import ClassifierMixin
from ..svm.base import BaseLibLinear
from ..svm import _liblinear

class LogisticRegression(BaseLibLinear, ClassifierMixin):
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
    that happens, try with a smaller eps parameter.

    References
    ----------
    LIBLINEAR -- A Library for Large Linear Classification
    http://www.csie.ntu.edu.tw/~cjlin/liblinear/
    """

    def __init__(self, penalty='l2', dual=False, eps=1e-4, C=1.0,
                 fit_intercept=True):

        super(LogisticRegression, self).__init__ (penalty=penalty,
            dual=dual, loss='lr', eps=eps, C=C,
            fit_intercept=fit_intercept)

    def predict_proba(self, T):
        """
        Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.
        """
        T = np.asanyarray(T, dtype=np.float64, order='C')
        probas = _liblinear.predict_prob_wrap(T, self.raw_coef_,
                                      self._get_solver_type(),
                                      self.eps, self.C,
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
