__all__ = [
    'ClassificationResults',
    'RegressionResults',
    'OneClassResults'
    ]

import libsvm
import utils

import numpy as N
from ctypes import *

class Results:
    def __init__(self, dtype, model):
        self.dtype = dtype
        self.model = model
        model = model.contents
        self.svm_type = model.param.svm_type
        self.nr_class = model.nr_class

    def __del__(self):
        libsvm.svm_destroy_model(self.model)
