__all__ = [
    'LibSvmOneClassModel'
    ]

from ctypes import cast, POINTER, byref, c_double

from model import LibSvmModel
import libsvm

class LibSvmOneClassResults:
    def __init__(self, model, dataset):
        self.model = model
        self.dataset = dataset
        model = model.contents
        self.rho = model.rho[0]
        self.sv_coef = model.sv_coef[0][:model.l]

    def __del__(self):
        libsvm.svm_destroy_model(self.model)

    def predict(self, dataset):
        """
        This function returns a list of boolean values which indicate
        whether the test vectors form part of the distribution.
        """
        def p(x):
            xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
            return libsvm.svm_predict(self.model, xptr) > 0
        return map(p, dataset.data)

    def predict_values(self, dataset):
        """
        This function returns a list of floating point values which
        indicate whether the test vectors form part of the
        distribution.

        A positive value indicates that the test vector is part of the
        distribution, while a non-positive value indicates that is is
        not.
        """
        def p(x):
            xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
            v = c_double()
            libsvm.svm_predict_values(self.model, xptr, byref(v))
            return v.value
        return map(p, dataset.data)

class LibSvmOneClassModel(LibSvmModel):
    """
    A model for distribution estimation (one-class SVM).

    See also: Scholkopf, et al. Estimating the Support of a
    High-Dimensional Distribution.
    """

    Results = LibSvmOneClassResults

    def __init__(self, kernel, nu=0.5, **kwargs):
        """
        Parameters:

        - `nu`: XXX
        """
        LibSvmModel.__init__(self, kernel, **kwargs)
        self.nu = nu
        self.param.svm_type = libsvm.ONE_CLASS
        self.param.nu = nu
