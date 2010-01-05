__all__ = [
    'LinearKernel',
    'PolynomialKernel',
    'RBFKernel',
    'SigmoidKernel',
    'CustomKernel'
    ]

def svm_node_dot(x, y):
    # associate node indexes with array indexes
    xidx = dict(zip(x['index'][:-1],range(0,len(x))))
    yidx = dict(zip(y['index'][:-1],range(0,len(y))))
    # indexes in either vector
    indexes = N.unique(N.hstack([x['index'],y['index']]))
    z = 0.
    for j in indexes:
        if j in xidx and j in yidx:
            # dot if index is present in both vectors
            z += x['value'][xidx[j]]*y['value'][yidx[j]]
    return z

class LinearKernel:
    def __init__(self, dot=svm_node_dot):
        self.dot = dot

    def __call__(self, x, y):
        return self.dot(x, y)

class PolynomialKernel:
    def __init__(self, degree, gamma, coef0, dot=svm_node_dot):
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
    def __init__(self, gamma, dot=svm_node_dot):
        self.gamma = gamma
        self.dot = dot

    def __call__(self, x, y):
        z = self.dot(x, x) + self.dot(y, y) - 2*self.dot(x, y)
        return N.exp(-self.gamma*z)

class SigmoidKernel:
    def __init__(self, gamma, coef0, dot=svm_node_dot):
        self.gamma = gamma
        self.coef0 = coef0
        self.dot = dot

    def kernel_sigmoid(x, y, gamma, coef0):
        return N.tanh(self.gamma*self.dot(x, y)+self.coef0)

class CustomKernel:
    """
    XXX example CustomKernel(lambda x, y, d: d(x,y))
    """

    def __init__(self, f, dot=svm_node_dot):
        self.f = f
        self.dot = dot

    def __call__(self, x, y):
        return self.f(x, y, dot)
