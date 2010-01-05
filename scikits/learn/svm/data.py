__all__ = [
    'LinearData',
    'LinearOneClassData',
    ]

import libsvm

import numpy as N

class LinearData:
    def __init__(self):
        self.kernel_type = libsvm.LINEAR

    def convert_train_data(self, data):
        svm_data = []
        for label, x in data:
            y = N.empty((len(x)+1,), dtype=libsvm.svm_node_dtype)
            y['index'][:-1] = N.arange(1, len(x)+1)
            y['value'][:-1] = x
            y[-1] = (-1, 0.)
            svm_data.append((label, y))
        return svm_data

    def convert_test_data(self, x):
        y = N.empty((len(x)+1,), dtype=libsvm.svm_node_dtype)
        y['index'][:-1] = N.arange(1, len(x)+1)
        y['value'][:-1] = x
        y[-1] = (-1, 0.)
        return y

class LinearOneClassData:
    def __init__(self):
        self.kernel_type = libsvm.LINEAR

    def convert_train_data(self, data):
        svm_data = []
        for x in data:
            y = N.empty((len(x)+1,), dtype=libsvm.svm_node_dtype)
            y['index'][:-1] = 0
            y['value'][:-1] = x
            y[-1] = (-1, 0.)
            svm_data.append((0, y))
        return svm_data

    def convert_test_data(self, x):
        y = N.empty((len(x)+1,), dtype=libsvm.svm_node_dtype)
        y['index'][:-1] = 0
        y['value'][:-1] = x
        y[-1] = (-1, 0.)
        return y
