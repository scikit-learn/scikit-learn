from ctypes import POINTER, c_double, addressof, byref
from itertools import izip
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
        if n == 1:
            return v[0]
        else:
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
        self.is_compact = False
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

    def _predict_values_sparse(self, x, n):
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
                    sum = 0.
                    for k in range(ci):
                        sum += coef1[si + k] * kvalue[si + k]
                    for k in range(cj):
                        sum += coef2[sj + k] * kvalue[sj + k]
                    dec_values[p] = sum - self.rho[p]
                    p += 1
            return dec_values
        else:
            z = -self.rho
            for sv_coef, sv in zip(self.sv_coef, self.support_vectors):
                z += sv_coef * self.kernel(x, sv, svm_node_dot)
            return z

    def _predict_values_compact(self, x, n):
        if self.svm_type in [libsvm.C_SVC, libsvm.NU_SVC]:
            for i, sv in enumerate(self.support_vectors):
                kvalue = N.empty((len(self.support_vectors),))
                kvalue[i] = self.kernel(x, sv, svm_node_dot)
            return kvalue - self.rho
        else:
            sv = self.support_vectors[0]
            return self.kernel(x, sv, svm_node_dot) - self.rho

    def predict_values(self, x, n):
        if self.is_compact:
            return self._predict_values_compact(x, n)
        else:
            return self._predict_values_sparse(x, n)

    def predict_probability(self, x, n):
        raise NotImplementedError

    def _compact_svs(self, svs, coefs):
        maxlen = 0
        for sv in svs:
            maxlen = N.maximum(maxlen, sv['index'].max())
        csv = N.zeros((maxlen + 1,), libsvm.svm_node_dtype)
        csv['index'][:-1] = N.arange(1, maxlen + 1)
        csv['index'][-1] = -1
        for coef, sv in izip(coefs, svs):
            idx = sv['index'][:-1] - 1
            csv['value'][idx] += coef*sv['value'][:-1]
        return csv

    def compact(self):
        if self.svm_type in [libsvm.C_SVC, libsvm.NU_SVC]:
            compact_support_vectors = []
            for i in range(self.nr_class):
                for j in range(i + 1, self.nr_class):
                    si, sj = self.start[i], self.start[j]
                    ci, cj = self.nSV[i], self.nSV[j]
                    svi = self.support_vectors[si:si + ci]
                    svj = self.support_vectors[sj:sj + cj]
                    coef1 = self.sv_coef[j - 1][si:si + ci]
                    coef2 = self.sv_coef[i][sj:sj + cj]
                    csv = self._compact_svs(svi + svj, coef1 + coef2)
                    compact_support_vectors.append(csv)
            self.support_vectors = compact_support_vectors
        else:
            csv = self._compact_svs(self.support_vectors, self.sv_coef)
            self.support_vectors = [csv]
        self.is_compact = True
