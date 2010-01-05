import numpy as N

from numpy.testing import *
set_local_path('../..')

from svm.dataset import convert_to_svm_node
from svm.dataset import *
from svm.libsvm import svm_node_dtype

restore_path()

class test_dataset(NumpyTestCase):
    def check_convert_dict(self):
        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node({}), x)

        x = N.array([(1,2.),(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node({1:2.}), x)

        x = N.array([(1,2.),(3,4.),(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node({1:2.,3:4.}), x)

        # check for positive indexes
        self.assertRaises(AssertionError, convert_to_svm_node, {0:0.})

    def check_convert_list(self):
        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node([]), x)

        x = N.array([(1,2.),(3,4.),(-1,0.)], dtype=svm_node_dtype)
        # check that indexes are sorted
        assert_array_equal(convert_to_svm_node([(3,4.),(1,2.)]), x)

        # check for unique indexes
        self.assertRaises(AssertionError,
                          convert_to_svm_node, [(1,0.),(1,0.)])

    def check_convert_array(self):
        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node(N.empty(0)), x)

        x = N.array([(1,1.),(2,2.),(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node(N.arange(1,3)), x)

    def check_regression(self):
        data = [(1.0, N.arange(5))]
        dataset = LibSvmRegressionDataSet(data)
        self.assertAlmostEqual(dataset.gamma, 0.2)

    def check_classification(self):
        data = [(1, N.arange(4)), (2, N.arange(10))]
        dataset = LibSvmClassificationDataSet(data)
        self.assertAlmostEqual(dataset.gamma, 0.1)
        self.assert_(1 in dataset.labels)
        self.assert_(2 in dataset.labels)

    def check_oneclass(self):
        data = [N.arange(2)]
        dataset = LibSvmOneClassDataSet(data)
        self.assertAlmostEqual(dataset.gamma, 0.5)

if __name__ == '__main__':
    NumpyTest().run()
