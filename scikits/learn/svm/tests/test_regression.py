import sys
import os
sys.path.insert(0, '..')
sys.path.insert(0, os.path.join('..','..'))

from numpy.testing import *
import numpy as N

from svm.regression import *
from svm.dataset import LibSvmRegressionDataSet
from svm.dataset import LibSvmTestDataSet
from svm.kernel import *

class test_regression(NumpyTestCase):
    def check_basics(self):
        Model = LibSvmEpsilonRegressionModel
        Kernel = LinearKernel()
        Model(Kernel)
        Model(Kernel, epsilon=0.1)
        Model(Kernel, cost=1.0)
        model = Model(Kernel, shrinking=False)
        self.assert_(not model.shrinking)

    def check_epsilon_train(self):
        y = [10., 20., 30., 40.]
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        traindata = LibSvmRegressionDataSet(zip(y, x))
        testdata = LibSvmTestDataSet(x)

        Model = LibSvmEpsilonRegressionModel
        model = Model(LinearKernel())

        results = model.fit(traindata)
        results.predict(testdata)
        results.get_svr_probability()

    def check_epsilon_more(self):
        y = [0.0, 1.0, 1.0, 2.0]
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        epsilon = 0.1
        cost = 10.0
        traindata = LibSvmRegressionDataSet(zip(y, x))
        testdata = LibSvmTestDataSet(x)

        kernels = [
            LinearKernel(),
            PolynomialKernel(3, traindata.gamma, 0.0),
            RBFKernel(traindata.gamma)
            ]
        expected_ys = [
            N.array([0.1, 1.0, 1.0, 1.9]),
            N.array([0.24611273, 0.899866638, 0.90006681, 1.90006681]),
            N.array([0.1, 1.0, 1.0, 1.9])
            ]

        for kernel, expected_y in zip(kernels, expected_ys):
            model = LibSvmEpsilonRegressionModel(kernel, epsilon, cost)
            results = model.fit(traindata)
            predictions = results.predict(testdata)
            assert_array_almost_equal(predictions, expected_y)

    def check_nu_train(self):
        pass

if __name__ == '__main__':
    NumpyTest().run()
