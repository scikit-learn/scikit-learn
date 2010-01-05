__all__ = [
    'LibSvmEpsilonRegressionModel',
    'LibSvmNuRegressionModel'
    ]

from model import LibSvmModel
import libsvm

"""
class RegressionResults(Results):
    def __init__(self, model):
        Results.__init__(self, model)
        model = model.contents
        self.rho = model.rho[0]
        self.sv_coef = model.sv_coef[0][:model.l]

    def predict(self, dataset):
        #x = self.dtype.convert_test_data(x)
        #xptr = utils.array_as_ctype(x, libsvm.svm_node)
        #return libsvm.svm_predict(self.model, xptr)
        raise NotImplementedError
"""

class LibSvmEpsilonRegressionModel(LibSvmModel):
    def __init__(self, kernel, epsilon=0.1, cost=1.0, **kwargs):
        LibSvmModel.__init__(self, libsvm.EPSILON_SVR, kernel, **kwargs)
        self.epsilon = epsilon
        self.cost = cost
        self.param.p = epsilon
        self.param.C = cost
        self.param.probability = 1

class LibSvmNuRegressionModel(LibSvmModel):
    def __init__(self, kernel, nu=0.5, cost=1.0, **kwargs):
        LibSvmModel.__init__(self, libsvm.NU_SVR, kernel, **kwargs)
        self.nu = nu
        self.cost = cost
        self.param.nu = nu
        self.param.C = cost
        self.param.probability = 1
