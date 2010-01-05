import numpy as N

from numpy.testing import *
set_local_path('../..')

import svm.libsvm as libsvm

restore_path()

class test_libsvm(NumpyTestCase):
    def check_svm_node(self):
        node = libsvm.svm_node()
        node = N.empty((), dtype=libsvm.svm_node_dtype)
        node = N.empty((1,), dtype=libsvm.svm_node_dtype)
        node[0]['index'] = 123
        node[0]['value'] = 456.
        assert_equal(node[0][0], 123)
        assert_equal(node[0][1], 456.)

    def check_svm_parameter(self):
        param = libsvm.svm_parameter()
        param.degree = 3
        param.gamma = 1.0

    def check_svm_problem(self):
        problem = libsvm.svm_problem()
        problem.l = 123

    def check_svm_model(self):
        model = libsvm.svm_model()
        model.nr_class = 123
        param = libsvm.svm_parameter()
        param.degree = 3
        model.param = param
        assert_equal(model.param.degree, 3)

if __name__ == '__main__':
    NumpyTest().run()
