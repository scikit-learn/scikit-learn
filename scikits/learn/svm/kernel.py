__all__ = [
    'LinearKernel',
    'PolynomialKernel',
    'RBFKernel',
    'SigmoidKernel',
    'CustomKernel'
    ]

class LinearKernel:
    def __init__(self, dot):
        self.dot = dot

    def __call__(self, x, y):
        return self.dot(x, y)

class PolynomialKernel:
    def __init__(self, degree, gamma, coef0, dot):
        self.degree = degree
        self.gamma = gamma
        self.coef0 = coef0
        self.dot = dot

    def __call__(self, x, y):
        base = self.gamma*self.dot(x, y) + self.coef0
        tmp = base
        ret = 1.
        t = self.degree
        while t > 0:
            if t % 2 == 1: ret *= tmp
            tmp *= tmp
            t /= 2
        return ret

class RBFKernel:
    def __init__(self, gamma, dot):
        self.gamma = gamma
        self.dot = dot

    def __call__(self, x, y):
        z = self.dot(x, x) + self.dot(y, y) - 2*self.dot(x, y)
        return N.exp(-self.gamma*z)

class SigmoidKernel:
    def __init__(self, gamma, coef0, dot):
        self.gamma = gamma
        self.coef0 = coef0
        self.dot = dot

    def kernel_sigmoid(x, y, gamma, coef0):
        return N.tanh(self.gamma*self.dot(x, y)+self.coef0)

class CustomKernel:
    def __init__(self, f, dot):
        self.f = f
        self.dot = dot

    def __call__(self, x, y):
        return self.f(x, y, dot)
