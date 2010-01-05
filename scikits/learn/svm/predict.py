from ctypes import cast, POINTER, c_double
import numpy as N

import libsvm

__all__ = [
    'LibSvmPredictor',
    'LibSvmPrecomputedPredictor',
    'LibSvmSparsePredictor',
    'PrecomputedSparsePredictor',
    'DensePredictor'
    ]

class LibSvmPredictor:
    def __init__(self, model, dataset):
        self.model = model

    def __del__(self):
        libsvm.svm_destroy_model(self.model)

    def predict(self, x):
        xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
        return libsvm.svm_predict(self.model, xptr)

    def predict_values(self, x, n):
        xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
        v = N.empty((n,), dtype=N.float64)
        vptr = cast(v.ctypes.data, POINTER(c_double))
        libsvm.svm_predict_values(self.model, xptr, vptr)
        return v

    def predict_probability(self, x, n):
        xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
        pe = N.empty((n,), dtype=N.float64)
        peptr = cast(pe.ctypes.data, POINTER(c_double))
        label = libsvm.svm_predict_probability(self.model, xptr, peptr)
        return label, pe

class LibSvmPrecomputedPredictor:
    def __init__(self, model, dataset):
        raise NotImplementedError

    def predict(self):
        raise NotImplementedError

    def predict_values(self):
        raise NotImplementedError

    def predict_probability(self):
        raise NotImplementedError

class LibSvmSparsePredictor:
    def __init__(self, model, dataset):
        raise NotImplementedError

    def predict_values(self):
        raise NotImplementedError

    def predict_probability(self):
        raise NotImplementedError

class PrecomputedSparsePredictor:
    def __init__(self, model, dataset):
        raise NotImplementedError

    def predict_values(self):
        raise NotImplementedError

    def predict_probability(self):
        raise NotImplementedError

class DensePredictor:
    def __init__(self, model, dataset):
        raise NotImplementedError

    def predict_values(self):
        raise NotImplementedError

    def predict_probability(self):
        raise NotImplementedError
