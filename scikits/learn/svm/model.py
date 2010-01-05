__all__ = [
    'Model'
    ]

import libsvm
import utils

import numpy as N
from ctypes import *

class Model:
    def __init__(self, dtype, shrinking=True, cache_size=40, tol=0.001):
        self.dtype = dtype
        self.shrinking = shrinking
        self.cache_size = cache_size
        self.tol = tol

    def fit(self, data):
        svm_data = self.dtype.convert_train_data(data)
        # libsvm requires data to be sorted by label
        svm_data.sort(cmp=lambda x, y: cmp(x[0], y[0]))
        param = self.setup_svm_parameter(svm_data)

        # XXX find better way to keep x and y references
        problem, x, y = self.setup_svm_problem(svm_data)

        self.check_problem_param(problem, param)
        model = libsvm.svm_train(problem, param)
        self.results = self.Results(self.dtype, model)

        # XXX find better way to keep svm_data reference
        self.results.svm_data = svm_data

        return self.results

    def predict(self, x):
        return self.results.predict(svm_data)

    def setup_svm_parameter(self, svm_data):
        param = libsvm.svm_parameter()
        param.svm_type = getattr(self, 'svm_type')
        param.kernel_type = getattr(self.dtype, 'kernel_type')
        param.degree = getattr(self.dtype, 'degree', 0)
        if hasattr(self.dtype, 'gamma') and self.dtype.gamma is None:
            maxlen = 0
            for x in svm_data:
                maxlen = max(maxlen, x[1]['index'][:-1].max())
            param.gamma = 1.0/maxlen
        else:
            param.gamma = getattr(self.dtype, 'gamma', 0.0)
        param.coef0 = getattr(self.dtype, 'coef0', 0)
        param.cache_size = getattr(self, 'cache_size')
        param.eps = getattr(self, 'tol')
        param.C = getattr(self, 'cost', 0.0)
        # XXX nr_weight, weight_label, weight
        param.nr_weight = 0
        # XXX setting these to None zeros svm_type
        ###param.weight_label = None
        ###param.weight = None
        param.nu = getattr(self, 'nu', 0.0)
        param.p = getattr(self, 'epsilon', 0.0)
        param.shrinking = getattr(self, 'shrinking')
        param.probability = 0
        return param

    def setup_svm_problem(self, svm_data):
        problem = libsvm.svm_problem()
        problem.l = len(svm_data)
        y = (c_double*problem.l)()
        x = (POINTER(libsvm.svm_node)*problem.l)()
        for i, (label, node) in enumerate(svm_data):
            y[i] = label
            x[i] = cast(utils.addressof_array(node),
                        POINTER(libsvm.svm_node))
        problem.x = cast(addressof(x), POINTER(POINTER(libsvm.svm_node)))
        problem.y = cast(addressof(y), POINTER(c_double))
        return problem, x, y

    def check_problem_param(self, problem, param):
        error_msg = libsvm.svm_check_parameter(problem, param)
        if error_msg:
            raise ValueError, error_msg
