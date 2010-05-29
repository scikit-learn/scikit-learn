import numpy as np
from . import liblinear


class LogisticRegression(object):
    """
    Logistic Regression.

    Implements L1 and L2 regularized logistic regression.

    Parameters
    ----------
    X : array-like, shape = [nsamples, nfeatures]
        Training vector, where nsamples in the number of samples and
        nfeatures is the number of features.
    Y : array, shape = [nsamples]
        Target vector relative to X

    penalty : string, 'l1' or 'l2'
        Used to specify the norm used in the penalization

    C : float
        Specifies the strength of the regularization. The smaller it is
        the bigger in the regularization.

    intercept : bool, default: True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added the decision function

    Attributes
    ----------

    `coef_` : array, shape = [nclasses-1, nfeatures]
        Coefficient of the features in the decision function.

    `intercept_` : array, shape = [nclasses-1]
        intercept (a.k.a. bias) added to the decision function.
        It is available only when parametr intercept is set to True

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
    def __init__(self, penalty='l2', eps=1e-4, C=1.0, intercept=True):
        self.solver_type = self._penalties[penalty.lower()]
        self.eps = eps
        self.C = C
        if intercept:
            self.bias_ = 1.0
        else:
            self.bias_ = -1.0

    _penalties = {'l2': 0, 'l1' : 6}
    _weight_label = np.empty(0, dtype=np.int32)
    _weight = np.empty(0, dtype=np.float64)

    def fit(self, X, Y):
        X = np.asanyarray(X, dtype=np.float64, order='C')
        Y = np.asanyarray(Y, dtype=np.int32, order='C')
        self.raw_coef_, self.label_, self.bias_ = liblinear.train_wrap(X,
                                          Y, self.solver_type, self.eps, self.bias_,
                                          self.C, 0,
                                          self._weight_label,
                                          self._weight)
        return self

    def predict(self, T):
        T = np.asanyarray(T, dtype=np.float64, order='C')
        return liblinear.predict_wrap(T, self.raw_coef_, self.solver_type,
                                      self.eps, self.C,
                                      self._weight_label,
                                      self._weight, self.label_,
                                      self.bias_)

    def predict_proba(self, T):
        T = np.asanyarray(T, dtype=np.float64, order='C')
        return liblinear.predict_prob_wrap(T, self.raw_coef_, self.solver_type,
                                      self.eps, self.C,
                                      self._weight_label,
                                      self._weight, self.label_,
                                      self.bias_)

    @property
    def intercept_(self):
        if self.bias_ > 0:
            return self.raw_coef_[:,-1]
        return 0.0
            

    @property
    def coef_(self):
        if self.bias_ > 0:
            return self.raw_coef_[:,:-1]
        return self.raw_coef_

