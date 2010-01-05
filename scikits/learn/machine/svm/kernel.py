import numpy as N

import libsvm

__all__ = [
    'LinearKernel',
    'PolynomialKernel',
    'RBFKernel',
    'SigmoidKernel',
    'CustomKernel'
    ]

class LinearKernel:
    """Linear kernel: u'*v."""
    def __init__(self):
        self.kernel_type = libsvm.LINEAR

    def __call__(self, x, y):
        x = N.atleast_2d(x)
        y = N.atleast_2d(y)
        return N.dot(x, y.T)

    def compact(self, *args):
        return self

class PolynomialKernel:
    """Polynomial kernel: (gamma*u'*v + coef0)^degree."""

    def __init__(self, degree, gamma, coef0):
        self.kernel_type = libsvm.POLY
        self.degree = degree
        self.gamma = gamma
        self.coef0 = coef0

    def __call__(self, x, y):
        x = N.atleast_2d(x)
        y = N.atleast_2d(y)
        base = self.gamma * N.dot(x, y.T) + self.coef0
        tmp = base
        ret = 1.0
        t = self.degree
        while t > 0:
            if t % 2 == 1: ret *= tmp
            tmp *= tmp
            t /= 2
        return ret

    def __repr__(self):
        return '<PolynomialKernel: degree=%d, gamma=%.4f, coef0=%.4f>' % \
            (self.degree, self.gamma, self.coef0)

    def compact(self, *args):
        raise NotImplementedError, \
            'model compaction for PolynomialKernel not implemented'

class RBFKernel:
    """Radial basis function kernel: exp(-gamma*|u-v|^2)"""

    def __init__(self, gamma):
        self.kernel_type = libsvm.RBF
        self.gamma = gamma
        self.__call__ = self.evaluate

    def evaluate(self, x, y):
        x = N.atleast_2d(x)
        y = N.atleast_2d(y)
        xnorm = N.atleast_2d(N.sum(x*x, axis=1))
        ynorm = N.atleast_2d(N.sum(y*y, axis=1))
        z = xnorm + ynorm - 2 * N.atleast_2d(N.dot(x, y.T).squeeze())
        return N.exp(-self.gamma * z)

    def evaluate_compact(self, x, y):
        raise NotImplementedError

    def __repr__(self):
        return '<RBFKernel: gamma=%.4f>' % (self.gamma,)

    def compact(self, *args):
        raise NotImplementedError, \
            'model compaction for RBFKernel not implemented'

class SigmoidKernel:
    """Sigmoid kernel: tanh(gamma*u'*v + coef0)"""

    def __init__(self, gamma, coef0):
        self.kernel_type = libsvm.SIGMOID
        self.gamma = gamma
        self.coef0 = coef0

    def __call__(self, x, y):
        x = N.atleast_2d(x)
        y = N.atleast_2d(y)
        return N.tanh(self.gamma * N.dot(x, y.T) + self.coef0)

    def __repr__(self):
        return '<SigmoidKernel: gamma=%.4f, coef0=%.4f>' % \
            (self.gamma, self.coef0)

    def compact(self, *args):
        raise NotImplementedError, \
            'model compaction for SigmoidKernel not implemented'

class CustomKernel:
    """Custom kernel: any callable"""

    def __init__(self, f):
        self.kernel_type = libsvm.PRECOMPUTED
        self.f = f

    def __call__(self, x, y):
        x = N.atleast_2d(x)
        y = N.atleast_2d(y)
        return self.f(x, y)

    def __repr__(self):
        return '<CustomKernel: %s>' % str(self.f)

    def compact(self, *args):
        raise NotImplementedError, \
            'model compaction for CustomKernel not implemented'
