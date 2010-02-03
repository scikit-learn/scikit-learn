from unittest import TestCase
from numpy.testing import assert_array_equal
import numpy as N

from ..dataset import OneClassDataSet, TestDataSet
from ..kernel import Linear, Polynomial, RBF
from ..oneclass import OneClassModel
from ..predict import PythonPredictor

class TestOneClass(TestCase):
    def test_basics(self):
        ModelType = OneClassModel
        kernel = Linear()
        ModelType(kernel)
        ModelType(kernel, nu=1.0)

    def _make_basic_datasets(self):
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        traindata = OneClassDataSet(x)
        testdata = TestDataSet(x)
        return traindata, testdata

    def test_train(self):
        traindata, testdata = self._make_basic_datasets()
        model = OneClassModel(Linear())
        results = model.fit(traindata)
        p = results.predict(testdata)
        assert_array_equal(p, [False, False, False, True])
        v = results.predict_values(testdata)
        assert_array_equal(v, [-0.5, 0.0, 0.0, 0.5])

    def test_more(self):
        traindata, testdata = self._make_basic_datasets()
        ModelType = OneClassModel
        nu = 0.5
        kernels = [
            Linear(),
            Polynomial(3, traindata.gamma, 0.0),
            RBF(traindata.gamma)
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

    def test_compact(self):
        traindata, testdata = self._make_basic_datasets()
        model = OneClassModel(Linear())
        results = model.fit(traindata, PythonPredictor)
        refv = results.predict_values(testdata)
        results.compact()
        v = results.predict_values(testdata)
        assert_array_equal(refv, v)
