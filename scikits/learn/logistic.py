import numpy as np
from . import liblinear

class Logistic(object):

    def __init__(self, eps=1e-4, C=1.0):
        self.eps = eps
        self.C = C

    _weight_label = np.empty(0, dtype=np.int)
    _weight = np.empty(0, dtype=np.float64)

    def fit(self, X, Y):
        X = np.asanyarray(X, dtype=np.float64, order='C')
        Y = np.asanyarray(Y, dtype=np.int, order='C')
        self.coef_, self.label_, self.bias_ = liblinear.train_wrap(X,
                                          Y, 0, self.eps, 1.0,
                                          self.C, 0,
                                          self._weight_label,
                                          self._weight)

    def predict(self, T):
        T = np.asanyarray(T, dtype=np.float64, order='C')
        return liblinear.predict_wrap(T, self.coef_, 0, self.eps,
                                      self.C, self._weight_label,
                                      self._weight, self.label_,
                                      self.bias_)
