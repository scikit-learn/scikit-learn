from numpy.testing import *
import numpy as N

from .. import libsvm as libsvm


class TestLibSvm(TestCase):
    def test_svm_node(self):
        node = libsvm.svm_node()
        node = N.empty((), dtype=libsvm.svm_node_dtype)
        node = N.empty((1,), dtype=libsvm.svm_node_dtype)
        node[0]['index'] = 123
        node[0]['value'] = 456.
        assert_equal(node[0][0], 123)
        assert_equal(node[0][1], 456.)

    def test_svm_parameter(self):
        param = libsvm.svm_parameter()
        param.degree = 3
        param.gamma = 1.0

    def test_svm_problem(self):
        problem = libsvm.svm_problem()
        problem.l = 123

    def test_svm_model(self):
        model = libsvm.svm_model()
        model.nr_class = 123
        param = libsvm.svm_parameter()
        param.degree = 3
        model.param = param
        assert_equal(model.param.degree, 3)

if __name__ == '__main__':
    run_module_suite()
