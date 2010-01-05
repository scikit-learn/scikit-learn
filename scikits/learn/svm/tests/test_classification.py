from numpy.testing import *
import numpy as N

from svm.classification import *
from svm.dataset import LibSvmClassificationDataSet
from svm.dataset import LibSvmTestDataSet
from svm.kernel import *

class test_classification(NumpyTestCase):
    def check_basics(self):
        Model = LibSvmCClassificationModel
        Kernel = LinearKernel()
        Model(Kernel)
        Model(Kernel, cost=1.0)
        weights = [(2, 10.0), (1, 20.0), (0, 30.0)]
        Model(Kernel, weights=weights)
        Model(Kernel, 1.0, weights)
        model = Model(Kernel, cost=1.0, weights=weights)

    def check_c_train(self):
        labels = [0, 1, 1, 2]
        x = [N.array([0, 0]),
             N.array([0, 1]),
             N.array([1, 0]),
             N.array([1, 1])]
        traindata = LibSvmClassificationDataSet(zip(labels, x))

        Model = LibSvmCClassificationModel
        model = Model(RBFKernel(traindata.gamma))
        results = model.fit(traindata)

        testdata = LibSvmTestDataSet(x)
        results.predict(testdata)
        results.predict_values(testdata)

    def check_nu_train(self):
        pass

if __name__ == '__main__':
    NumpyTest().run()
