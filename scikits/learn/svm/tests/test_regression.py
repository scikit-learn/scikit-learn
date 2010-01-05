from numpy.testing import *
import numpy as N

from svm.regression import *
from svm.dataset import LibSvmRegressionDataSet
from svm.kernel import LinearKernel

class test_regression(NumpyTestCase):
    def check_basics(self):
        Model = LibSvmEpsilonRegressionModel
        Kernel = LinearKernel()
        Model(Kernel)
        Model(Kernel, epsilon=0.1)
        Model(Kernel, cost=1.0)
        model = Model(Kernel, shrinking=False)
        self.assert_(not model.shrinking)

    def check_epsilon(self):
        y = [10., 20., 30., 40.]
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        dataset = LibSvmRegressionDataSet(zip(y, x))

        Model = LibSvmEpsilonRegressionModel
        model = Model(LinearKernel())
        model.fit(dataset)

if __name__ == '__main__':
    NumpyTest().run()
