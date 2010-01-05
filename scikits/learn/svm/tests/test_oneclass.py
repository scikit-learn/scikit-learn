from numpy.testing import *
import numpy as N

from svm.oneclass import *
from svm.dataset import LibSvmOneClassDataSet
from svm.dataset import LibSvmTestDataSet
from svm.kernel import *

class test_oneclass(NumpyTestCase):
    def check_basics(self):
        Model = LibSvmOneClassModel
        Kernel = LinearKernel()
        Model(Kernel)
        Model(Kernel, nu=1.0)

    def check_train(self):
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        traindata = LibSvmOneClassDataSet(x)

        Model = LibSvmOneClassModel
        model = Model(LinearKernel())
        results = model.fit(traindata)

        testdata = LibSvmTestDataSet(x)
        results.predict(testdata)
        results.predict_values(testdata)

    def check_more(self):
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        traindata = LibSvmOneClassDataSet(x)
        nu = 0.5
        testdata = LibSvmTestDataSet(x)

        kernels = [
            LinearKernel(),
            PolynomialKernel(3, traindata.gamma, 0.0),
            RBFKernel(traindata.gamma)
            ]
        expected_preds = [
            [False, False, False, True],
            [False, False, False, True],
            [True, False, False, False]
            ]

        for kernel, expected_pred in zip(kernels, expected_preds):
            model = LibSvmOneClassModel(kernel, nu)
            results = model.fit(traindata)
            pred = results.predict(testdata)
            self.assertEqual(results.predict(testdata), expected_pred)
            values = results.predict_values(testdata)
            for p, v in zip(pred, values):
                self.assertEqual(v > 0, p)

if __name__ == '__main__':
    NumpyTest().run()
