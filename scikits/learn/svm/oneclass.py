from model import Model
from results import Results
import libsvm
import utils

from ctypes import c_double, byref

class OneClassResults(Results):
    def __init__(self, dtype, model):
        Results.__init__(self, dtype, model)
        model = model.contents
        self.rho = model.rho[0]
        self.sv_coef = model.sv_coef[0][:model.l]

    def predict_values(self, x):
        x = self.dtype.convert_test_data(x)
        v = c_double(0.0)
        xptr = utils.array_as_ctype(x, libsvm.svm_node)
        libsvm.svm_predict_values(self.model, xptr, byref(v))
        return v.value

    def predict(self, x):
        x = self.dtype.convert_test_data(x)
        xptr = utils.array_as_ctype(x, libsvm.svm_node)
        v = libsvm.svm_predict(self.model, xptr)
        return int(v)

class OneClassModel(Model):
    """
    A model for distribution estimation (one-class SVM).

    See also: Scholkopf, et al.: Estimating the Support of a
    High-Dimensional Distribution
    """

    Results = OneClassResults

    def __init__(self, dtype, nu=0.5, **kwargs):
        Model.__init__(self, dtype, **kwargs)
        self.svm_type = libsvm.ONE_CLASS
        self.nu = nu
