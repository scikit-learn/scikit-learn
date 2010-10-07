import numpy as np

from .base import ClassifierMixin
from .svm.base import BaseLibLinear
from .svm import _liblinear

class LogisticRegression(BaseLibLinear, ClassifierMixin):
    """
    Logistic Regression.

    Implements L1 and L2 regularized logistic regression.

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.
    Y : array, shape = [n_samples]
        Target vector relative to X

    penalty : string, 'l1' or 'l2'
        Used to specify the norm used in the penalization

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

    Methods
    -------
    fit(X, Y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    See also
    --------
    LinearSVC

    References
    ----------
    LIBLINEAR -- A Library for Large Linear Classification
    http://www.csie.ntu.edu.tw/~cjlin/liblinear/
    """

    def __init__(self, penalty='l2', eps=1e-4, C=1.0, fit_intercept=True):
        super(LogisticRegression, self).__init__ (penalty=penalty, loss='lr',
            dual=False, eps=eps, C=C, fit_intercept=fit_intercept)

    def predict_proba(self, T):
        T = np.asanyarray(T, dtype=np.float64, order='C')
        return _liblinear.predict_prob_wrap(T, self.raw_coef_, self._get_solver_type(),
                                      self.eps, self.C,
                                      self._weight_label,
                                      self._weight, self.label_,
                                      self._get_bias())
