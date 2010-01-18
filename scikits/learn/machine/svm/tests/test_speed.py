from unittest import TestCase
import numpy as N

from ..classification import ClassificationModel, CClassificationModel
from ..dataset import TestDataSet, ClassificationDataSet
from ..kernel import Linear
from ..predict import PythonPredictor


class TestClassificationSpeed(TestCase):
    def check_large_test_dataset(self):
        x = N.random.randn(150, 3)

        # XXX shows bug where we can't get any support vectors 
        #x = N.random.randn(4, 2)

        #x = N.random.randn(10, 3)

        labels = N.random.random_integers(1, 5, x.shape[0])
        #labels = N.random.random_integers(1, 2, x.shape[0])
        traindata = ClassificationDataSet(labels, x)
        #kernel = RBFKernel(traindata.gamma)
        kernel = Linear()
        #kernel = PolynomialKernel(2, 5, 10)
        model = CClassificationModel(kernel)
        #xdim, ydim, zdim = 1, 1, x.shape[1]
        xdim, ydim, zdim = 2, 2, x.shape[1]
        img = N.random.randn(xdim, ydim, zdim)
        testdata1 = TestDataSet(img.reshape(xdim*ydim, zdim))
        testdata2 = TestDataSet(list(img.reshape(xdim*ydim, zdim)))

        refresults = model.fit(traindata)
        refv1 = refresults.predict_values(testdata1)
        refv2 = refresults.predict_values(testdata2)

        results = model.fit(traindata, PythonPredictor)
        #v11 = results.predict_values(testdata1)
        #v12 = results.predict_values(testdata2)

        results.compact()
        v21 = results.predict_values(testdata1)
        #v22 = results.predict_values(testdata2)

        print refv1
        print refv2
        #print v11
        #print v12
        print v21
        #print v22
