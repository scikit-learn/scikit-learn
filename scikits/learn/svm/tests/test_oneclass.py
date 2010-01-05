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

    def _make_basic_datasets(self):
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        traindata = LibSvmOneClassDataSet(x)
        testdata = LibSvmTestDataSet(x)
        return traindata, testdata

    def check_train(self):
        traindata, testdata = self._make_basic_datasets()
        model = LibSvmOneClassModel(LinearKernel())
        results = model.fit(traindata)
        p = results.predict(testdata)
        assert_array_equal(p, [False, False, False, True])
        v = results.predict_values(testdata)
        assert_array_equal(v, [-0.5, 0.0, 0.0, 0.5])

    def check_more(self):
        traindata, testdata = self._make_basic_datasets()
        ModelType = LibSvmOneClassModel
        nu = 0.5
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
            results = model.fit(traindata)
            pred = results.predict(testdata)
            self.assertEqual(results.predict(testdata), expected_pred)
            values = results.predict_values(testdata)
            for p, v in zip(pred, values):
                self.assertEqual(v > 0, p)

    def check_compact(self):
        traindata, testdata = self._make_basic_datasets()
        model = LibSvmOneClassModel(LinearKernel())
        results = model.fit(traindata, LibSvmPythonPredictor)
        refv = results.predict_values(testdata)
        results.compact()
        v = results.predict_values(testdata)
        assert_array_equal(refv, v)

if __name__ == '__main__':
    NumpyTest().run()
