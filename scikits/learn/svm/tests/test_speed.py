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
        x = N.random.randn(150, 5)
        labels = N.random.random_integers(1, 5, x.shape[0])
        traindata = LibSvmClassificationDataSet(labels, x)

        kernel = RBFKernel(traindata.gamma)
        model = LibSvmCClassificationModel(kernel)
        results = model.fit(traindata, LibSvmPythonPredictor)
        results.compact()

        xdim, ydim = 32, 32
        img = N.random.randn(xdim, ydim, 3)
        testdata = LibSvmTestDataSet(img.reshape(xdim*ydim, 3))
        v = results.predict_values(testdata)

if __name__ == '__main__':
    NumpyTest().run()
