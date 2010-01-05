from numpy.testing import *
import numpy as N

set_local_path('../..')
from svm.dataset import *
from svm.kernel import *
from svm.dataset import convert_to_svm_node, svm_node_dot
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

class test_svm_node_dot(NumpyTestCase):
    def check_dot(self):
        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        self.assertAlmostEqual(svm_node_dot(x, x), 0.)

        x = N.array([(1,1.),(-1,0.)], dtype=svm_node_dtype)
        y = N.array([(2,2.),(-1,0.)], dtype=svm_node_dtype)
        self.assertAlmostEqual(svm_node_dot(x, y), 0.)

        x = N.array([(3,2.),(-1,0.)], dtype=svm_node_dtype)
        y = N.array([(3,2.),(-1,0.)], dtype=svm_node_dtype)
        self.assertAlmostEqual(svm_node_dot(x, y), 4.)

class test_precomputed_dataset(NumpyTestCase):
    def check_precompute(self):
        degree, gamma, coef0 = 4, 3.0, 2.0
        kernels = [
            LinearKernel(),
            PolynomialKernel(degree, gamma, coef0),
            RBFKernel(gamma),
            SigmoidKernel(gamma, coef0)
            ]
        y = N.random.randn(10)
        x = N.random.randn(len(y), 10)
        origdata = LibSvmRegressionDataSet(zip(y, x))

        for kernel in kernels:
            # calculate expected Gram matrix
            expt_grammat = N.empty((len(y),)*2)
            for i, xi in enumerate(x):
                for j, xj in enumerate(x):
                    expt_grammat[i, j] = kernel(xi, xj, N.dot)
            # get a new dataset containing the precomputed data
            pcdata = origdata.precompute(kernel)
            for i, row in enumerate(pcdata.grammat):
                valuerow = row[1:]['value']
                assert_array_almost_equal(valuerow, expt_grammat[i])

    def check_combine(self):
        kernel = LinearKernel()

        y1 = N.random.randn(10)
        x1 = N.random.randn(len(y1), 10)
        origdata = LibSvmRegressionDataSet(zip(y1, x1))
        pcdata = origdata.precompute(kernel)

        y2 = N.random.randn(5)
        x2 = N.random.randn(len(y2), x1.shape[1])
        moredata = LibSvmRegressionDataSet(zip(y2, x2))
        morepcdata = pcdata.combine(moredata)

        expt_grammat = N.empty((len(y1) + len(y2),)*2)
        x = N.vstack([x1,x2])
        for i, xi in enumerate(x):
            for j, xj in enumerate(x):
                expt_grammat[i, j] = kernel(xi, xj, N.dot)
        for i, row in enumerate(morepcdata.grammat):
            valuerow = row[1:]['value']
            assert_array_almost_equal(valuerow, expt_grammat[i])

if __name__ == '__main__':
    NumpyTest().run()
