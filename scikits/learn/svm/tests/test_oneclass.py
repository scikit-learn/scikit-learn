from numpy.testing import *
import numpy as N

set_local_path('../..')
from svm.dataset import LibSvmOneClassDataSet, LibSvmTestDataSet
from svm.kernel import *
from svm.oneclass import *
from svm.predict import *
restore_path()

class test_oneclass(NumpyTestCase):
    def check_basics(self):
        ModelType = LibSvmOneClassModel
        kernel = LinearKernel()
        ModelType(kernel)
        ModelType(kernel, nu=1.0)

    def check_train(self):
        ModelType = LibSvmOneClassModel
        ResultType = LibSvmOneClassResults
        PredictorType = LibSvmPredictor

        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        traindata = LibSvmOneClassDataSet(x)
        model = ModelType(LinearKernel())
        results = model.fit(traindata, ResultType, PredictorType)
        testdata = LibSvmTestDataSet(x)
        results.predict(testdata)
        results.predict_values(testdata)

    def check_more(self):
        ModelType = LibSvmOneClassModel
        ResultType = LibSvmOneClassResults
        PredictorType = LibSvmPredictor

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
            model = ModelType(kernel, nu)
            results = model.fit(traindata, ResultType, PredictorType)
            pred = results.predict(testdata)
            self.assertEqual(results.predict(testdata), expected_pred)
            values = results.predict_values(testdata)
            for p, v in zip(pred, values):
                self.assertEqual(v > 0, p)

if __name__ == '__main__':
    NumpyTest().run()
