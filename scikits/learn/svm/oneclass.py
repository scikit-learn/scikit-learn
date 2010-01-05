__all__ = [
    'LibSvmOneClassModel'
    ]

from ctypes import cast, POINTER

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
        def p(x):
            xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
            # XXX maybe we want to cast return value to int
            return libsvm.svm_predict(self.model, xptr)
        return map(p, dataset.data)

    # XXX predict_values might also be useful

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
