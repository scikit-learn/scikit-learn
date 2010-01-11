from unittest import TestCase
from numpy.testing import assert_array_equal, assert_array_almost_equal
import numpy as N

from ..dataset import OneClassDataSet, RegressionDataSet, \
     ClassificationDataSet
from ..kernel import Linear, Polynomial, RBF, Sigmoid
from ..dataset import convert_to_svm_node, svm_node_dot
from ..libsvm import svm_node_dtype


class TestDataset(TestCase):
    def test_convert_dict(self):
        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node({}), x)

        x = N.array([(1,2.),(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node({1:2.}), x)

        x = N.array([(1,2.),(3,4.),(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node({1:2.,3:4.}), x)

        # check for positive indexes
        self.assertRaises(AssertionError, convert_to_svm_node, {0:0.})

    def test_convert_list(self):
        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node([]), x)

        x = N.array([(1,2.),(3,4.),(-1,0.)], dtype=svm_node_dtype)
        # check that indexes are sorted
        assert_array_equal(convert_to_svm_node([(3,4.),(1,2.)]), x)

        # check for unique indexes
        self.assertRaises(AssertionError,
                          convert_to_svm_node, [(1,0.),(1,0.)])

    def test_convert_array(self):
        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node(N.empty(0)), x)

        x = N.array([(1,1.),(2,2.),(-1,0.)], dtype=svm_node_dtype)
        assert_array_equal(convert_to_svm_node(N.arange(1,3)), x)

    def test_regression(self):
        data = [(1.0, N.arange(5))]
        y = map(lambda x: x[0], data)
        x = map(lambda x: x[1], data)
        dataset = RegressionDataSet(y, x)
        self.assertAlmostEqual(dataset.gamma, 0.2)
        self.assertEqual(len(dataset), len(data))
        for i, x in enumerate(dataset):
            self.assertEqual(data[i][0], x[0])
            assert_array_equal(data[i][1], x[1]['value'][:-1])

    def test_classification(self):
        data = [(1, N.arange(4)), (2, N.arange(10))]
        labels = map(lambda x: x[0], data)
        x = map(lambda x: x[1], data)
        dataset = ClassificationDataSet(labels, x)
        self.assertAlmostEqual(dataset.gamma, 0.1)
        #self.assert_(1 in dataset.labels)
        #self.assert_(2 in dataset.labels)
        self.assertEqual(len(dataset), len(data))
        for i, x in enumerate(dataset):
            self.assertEqual(data[i][0], x[0])
            assert_array_equal(data[i][1], x[1]['value'][:-1])

    def test_oneclass(self):
        data = [N.arange(2)]
        dataset = OneClassDataSet(data)
        self.assertAlmostEqual(dataset.gamma, 0.5)
        self.assertEqual(len(dataset), len(data))
        for i, x in enumerate(dataset):
            assert_array_equal(data[i], x[1]['value'][:-1])

class TestSvmNodeDot(TestCase):
    def test_basics(self):
        kernel = Linear()

        x = N.array([(-1,0.)], dtype=svm_node_dtype)
        self.assertAlmostEqual(svm_node_dot(x, x, kernel), 0.)

        x = N.array([(1,1.),(-1,0.)], dtype=svm_node_dtype)
        y = N.array([(2,2.),(-1,0.)], dtype=svm_node_dtype)
        self.assertAlmostEqual(svm_node_dot(x, y, kernel), 0.)

        x = N.array([(3,2.),(-1,0.)], dtype=svm_node_dtype)
        y = N.array([(3,2.),(-1,0.)], dtype=svm_node_dtype)
        self.assertAlmostEqual(svm_node_dot(x, y, kernel), 4.)

class TestPrecomputedDataset(TestCase):
    def test_precompute(self):
        degree, gamma, coef0 = 4, 3.0, 2.0
        kernels = [
            Linear(),
            Polynomial(degree, gamma, coef0),
            RBF(gamma),
            Sigmoid(gamma, coef0)
            ]
        y = N.random.randn(10)
        x = N.random.randn(len(y), 10)
        origdata = RegressionDataSet(y, x)

        for kernel in kernels:
            # calculate expected Gram matrix
            expt_grammat = N.empty((len(y),)*2)
            for i, xi in enumerate(x):
                for j, xj in enumerate(x):
                    expt_grammat[i, j] = kernel(xi, xj)
            # get a new dataset containing the precomputed data
            pcdata = origdata.precompute(kernel)
            for i, row in enumerate(pcdata.grammat):
                valuerow = row[1:]['value']
                assert_array_almost_equal(valuerow, expt_grammat[i])

    def test_combine(self):
        kernel = Linear()

        y1 = N.random.randn(10)
        x1 = N.random.randn(len(y1), 10)
        origdata = RegressionDataSet(y1, x1)
        pcdata = origdata.precompute(kernel)

        y2 = N.random.randn(5)
        x2 = N.random.randn(len(y2), x1.shape[1])
        moredata = RegressionDataSet(y2, x2)
        morepcdata = pcdata.combine(moredata)

        expt_grammat = N.empty((len(y1) + len(y2),)*2)
        x = N.vstack([x1,x2])
        for i, xi in enumerate(x):
            for j, xj in enumerate(x):
                expt_grammat[i, j] = kernel(xi, xj)
        for i, row in enumerate(morepcdata.grammat):
            valuerow = row[1:]['value']
            assert_array_almost_equal(valuerow, expt_grammat[i])
