__all__ = [
    'LinearKernel',
    'PolynomialKernel',
    'RBFKernel',
    'SigmoidKernel',
    'CustomKernel'
    ]

import numpy as N

class LinearKernel:
    def __call__(self, x, y, dot):
        return dot(x, y)

class PolynomialKernel:
    def __init__(self, degree, gamma, coef0):
        self.degree = degree
        self.gamma = gamma
        self.coef0 = coef0

    def __call__(self, x, y, dot):
        base = self.gamma*dot(x, y) + self.coef0
        tmp = base
        ret = 1.0
        t = self.degree
        while t > 0:
            if t % 2 == 1: ret *= tmp
            tmp *= tmp
            t /= 2
        return ret

class RBFKernel:
    def __init__(self, gamma):
        self.gamma = gamma

    def __call__(self, x, y, dot):
        z = dot(x, x) + dot(y, y) - 2*dot(x, y)
        return N.exp(-self.gamma*z)

class SigmoidKernel:
    def __init__(self, gamma, coef0):
        self.gamma = gamma
        self.coef0 = coef0

    def __call__(self, x, y, dot):
        return N.tanh(self.gamma*dot(x, y)+self.coef0)

class CustomKernel:
    def __init__(self, f):
        self.f = f

    def __call__(self, x, y, dot):
        return self.f(x, y, dot)
