from ctypes import POINTER, c_double, addressof
import numpy as N

from dataset import svm_node_dot
import libsvm

__all__ = [
    'LibSvmPredictor',
    'LibSvmPythonPredictor'
    ]

class LibSvmPredictor:
    def __init__(self, model, dataset, kernel):
        self.model = model
        self.kernel = kernel
        modelc = model.contents
        if modelc.param.kernel_type == libsvm.PRECOMPUTED:
            self.dataset = dataset
            self.sv_ids = [int(modelc.SV[i][0].value)
                           for i in range(modelc.l)]
            self._transform_input = self._create_gramvec
        else:
            self._transform_input = lambda x: x

    def __del__(self):
        libsvm.svm_destroy_model(self.model)

    def _create_gramvec(self, x):
        gramvec = N.zeros((len(self.dataset)+1,),
                          dtype=libsvm.svm_node_dtype)
        for sv_id in self.sv_ids:
            sv = self.dataset[sv_id]
            gramvec[sv_id]['value'] = self.kernel(x, sv, svm_node_dot)
        return gramvec

    def predict(self, x):
        x = self._transform_input(x)
        xptr = x.ctypes.data_as(POINTER(libsvm.svm_node))
        return libsvm.svm_predict(self.model, xptr)

    def predict_values(self, x, n):
        x = self._transform_input(x)
        xptr = x.ctypes.data_as(POINTER(libsvm.svm_node))
        v = N.empty((n,), dtype=N.float64)
        vptr = v.ctypes.data_as(POINTER(c_double))
        libsvm.svm_predict_values(self.model, xptr, vptr)
        return v

    def predict_probability(self, x, n):
        if not self.model.contents.param.probability:
            raise ValueError, 'not a probability model'
        x = self._transform_input(x)
        xptr = x.ctypes.data_as(POINTER(libsvm.svm_node))
        pe = N.empty((n,), dtype=N.float64)
        peptr = pe.ctypes.data_as(POINTER(c_double))
        label = libsvm.svm_predict_probability(self.model, xptr, peptr)
        return label, pe

class LibSvmPythonPredictor:
    def __init__(self, model, dataset, kernel):
        self.kernel = kernel
        modelc = model.contents

        self.rho = modelc.rho[0]
        self.sv_coef = modelc.sv_coef[0][:modelc.l]
        self.svm_type = modelc.param.svm_type

        if modelc.param.kernel_type != libsvm.PRECOMPUTED:
            svptrs = [modelc.SV[i] for i in range(modelc.l)]
            support_vectors = [dataset.iddatamap[addressof(svptr[0])]
                               for svptr in svptrs]
        else:
            ids = [int(modelc.SV[i][0].value) for i in range(modelc.l)]
            support_vectors = [dataset[id] for id in ids]
        self.support_vectors = support_vectors

        libsvm.svm_destroy_model(model)

    def predict(self, x):
        if self.svm_type in [libsvm.C_SVC, libsvm.NU_SVC]:
            raise NotImplementedError
        else:
            return self.predict_values(x, 1)

    def predict_values(self, x, n):
        z = -self.rho
        # XXX possible optimization: izip
        for sv_coef, sv in zip(self.sv_coef, self.support_vectors):
            z += sv_coef * self.kernel(x, sv, svm_node_dot)
        return z

    def predict_probability(self, x, n):
        raise NotImplementedError
