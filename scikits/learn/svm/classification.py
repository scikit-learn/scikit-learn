from model import Model
from results import Results
import libsvm
import utils

import numpy as N
from ctypes import c_double

# XXX classification models need weights

class ClassificationResults(Results):
    def __init__(self, dtype, model):
        Results.__init__(self, dtype, model)
        model = model.contents
        self.labels = model.labels[:self.nr_class]
        self.rho = model.rho[:self.nr_class*(self.nr_class-1)/2]
        self.nSV = model.nSV[:self.nr_class]
        sv_coef = N.empty((self.nr_class-1, model.l), dtype=N.float64)
        for i, c in enumerate(model.sv_coef[:self.nr_class-1]):
            sv_coef[i,:] = c[:model.l]
        self.sv_coef = sv_coef

    def predict(self, x):
        x = self.dtype.convert_test_data(x)
        xptr = utils.array_as_ctype(x, libsvm.svm_node)
        return int(libsvm.svm_predict(self.model, xptr))

    def predict_values(self, x):
        x = self.dtype.convert_test_data(x)
        n = self.nr_class*(self.nr_class-1)/2
        v = N.empty((n,), dtype=N.float64)
        xptr = utils.array_as_ctype(x, libsvm.svm_node)
        vptr = utils.array_as_ctype(v, c_double)
        libsvm.svm_predict_values(self.model, xptr, vptr)
        count = 0
        d = {}
        for i in range(len(self.labels)):
            for j in range(i+1, len(self.labels)):
                d[self.labels[i], self.labels[j]] = v[count]
                d[self.labels[j], self.labels[i]] = -v[count]
                count += 1
        return d

class ClassificationModel(Model):
    Results = ClassificationResults

    def __init__(self, dtype, **kwargs):
        Model.__init__(self, dtype, **kwargs)

    def predict_values(self, x):
        return self.results.predict_values(svm_data)

class CSVCModel(ClassificationModel):
    def __init__(self, dtype, cost=1.0, **kwargs):
        ClassificationModel.__init__(self, dtype, **kwargs)
        self.svm_type = libsvm.C_SVC
        self.cost = cost

class NuSVCModel(ClassificationModel):
    def __init__(self, dtype, nu=0.5, **kwargs):
        ClassificationModel.__init__(self, dtype, **kwargs)
        self.svm_type = libsvm.NU_SVC
        self.nu = nu
