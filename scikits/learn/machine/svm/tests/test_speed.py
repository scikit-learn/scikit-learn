from numpy.testing import *
import numpy as N

set_local_path('../..')
from svm.classification import *
from svm.dataset import *
from svm.kernel import *
from svm.predict import *
restore_path()

class test_classification_speed(NumpyTestCase):
    def check_large_test_dataset(self):
        x = N.random.randn(150, 3)

        # XXX shows bug where we can't get any support vectors 
        #x = N.random.randn(4, 2)

        #x = N.random.randn(10, 3)

        labels = N.random.random_integers(1, 5, x.shape[0])
        #labels = N.random.random_integers(1, 2, x.shape[0])
        traindata = LibSvmClassificationDataSet(labels, x)
        #kernel = RBFKernel(traindata.gamma)
        kernel = LinearKernel()
        #kernel = PolynomialKernel(2, 5, 10)
        model = LibSvmCClassificationModel(kernel)
        #xdim, ydim, zdim = 1, 1, x.shape[1]
        xdim, ydim, zdim = 2, 2, x.shape[1]
        img = N.random.randn(xdim, ydim, zdim)
        testdata1 = LibSvmTestDataSet(img.reshape(xdim*ydim, zdim))
        testdata2 = LibSvmTestDataSet(list(img.reshape(xdim*ydim, zdim)))

        refresults = model.fit(traindata)
        refv1 = refresults.predict_values(testdata1)
        refv2 = refresults.predict_values(testdata2)

        results = model.fit(traindata, LibSvmPythonPredictor)
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

if __name__ == '__main__':
    NumpyTest().run()
