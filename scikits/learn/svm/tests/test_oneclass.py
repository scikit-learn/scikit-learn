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
        dataset = LibSvmOneClassDataSet(x)
        
        Model = LibSvmOneClassModel
        model = Model(LinearKernel())
        results = model.fit(dataset)

        testdata = LibSvmTestDataSet(x)
        results.predict(testdata)

if __name__ == '__main__':
    NumpyTest().run()
