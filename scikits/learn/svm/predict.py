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

    def compact(self):
        raise NotImplementedError

class LibSvmPythonPredictor:
    def __init__(self, model, dataset, kernel):
        self.kernel = kernel
        modelc = model.contents
        self.svm_type = modelc.param.svm_type
        if self.svm_type in [libsvm.C_SVC, libsvm.NU_SVC]:
            self.nr_class = modelc.nr_class
            self.labels = modelc.labels[:self.nr_class]
            nrho = self.nr_class * (self.nr_class - 1) / 2
            self.rho = modelc.rho[:nrho]
            self.sv_coef = [modelc.sv_coef[i][:modelc.l]
                            for i in range(self.nr_class - 1)]
            self.nSV = [modelc.nSV[i] for i in range(self.nr_class)]
            start = N.zeros((self.nr_class,), N.intc)
            for i in range(1, self.nr_class):
                start[i] = start[i - 1] + modelc.nSV[i - 1]
            self.start = start
        else:
            self.rho = modelc.rho[0]
            self.sv_coef = modelc.sv_coef[0][:modelc.l]

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
            n = self.nr_class * (self.nr_class - 1) / 2
            dec_values = self.predict_values(x, n)
            vote = N.zeros((self.nr_class,), N.intc)
            pos = 0
            for i in range(self.nr_class):
                for j in range(i + 1, self.nr_class):
                    if dec_values[pos] > 0:
                        vote[i] += 1
                    else:
                        vote[j] += 1
                    pos += 1
            return self.labels[vote.argmax()]
        else:
            return self.predict_values(x, 1)

    def predict_values(self, x, n):
        if self.svm_type in [libsvm.C_SVC, libsvm.NU_SVC]:
            kvalue = N.empty((len(self.support_vectors),))
            for i, sv in enumerate(self.support_vectors):
                kvalue[i] = self.kernel(x, sv, svm_node_dot)
            p = 0
            dec_values = N.empty((n,))
            for i in range(self.nr_class):
                for j in range(i + 1, self.nr_class):
                    sum = 0
                    si, sj = self.start[i], self.start[j]
                    ci, cj = self.nSV[i], self.nSV[j]
                    coef1 = self.sv_coef[j - 1]
                    coef2 = self.sv_coef[i]
                    sum = -self.rho[p]
                    for k in range(ci):
                        sum += coef1[si + k] * kvalue[si + k]
                    for k in range(cj):
                        sum += coef2[sj + k] * kvalue[sj + k]
                    dec_values[p] = sum
                    p += 1
            return dec_values
        else:
            z = -self.rho
            for sv_coef, sv in zip(self.sv_coef, self.support_vectors):
                z += sv_coef * self.kernel(x, sv, svm_node_dot)
            return z

    def predict_probability(self, x, n):
        raise NotImplementedError

    def compact(self):
        raise NotImplementedError
