from model import Model
from results import Results

import libsvm

class RegressionResults(Results):
    def __init__(self, dtype, model):
        Results.__init__(self, dtype, model)
        model = model.contents
        self.rho = model.rho[0]
        self.sv_coef = model.sv_coef[0][:model.l]

    def predict(self, x):
        x = self.dtype.convert_test_data(x)
        xptr = cast(utils.addressof_array(x), POINTER(libsvm.svm_node))
        return libsvm.svm_predict(self.model, xptr)

class EpsilonSVRModel(Model):
    Results = RegressionResults

    def __init__(self, dtype, cost=1.0, epsilon=0.1, **kwargs):
        Model.__init__(self, dtype, **kwargs)
        self.svm_type = libsvm.EPSILON_SVR
        self.cost = cost
        self.epsilon = epsilon

class NuSVRModel(Model):
    Results = RegressionResults

    def __init__(self, dtype, cost=1.0, nu=0.5, **kwargs):
        Model.__init__(self, dtype, **kwargs)
        self.svm_type = libsvm.NU_SVR
        self.cost = cost
        self.nu = nu
