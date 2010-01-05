from numpy.testing import *

# XXX remove this
import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))

import svm
import numpy as N

class test_oneclass(NumpyTestCase):
    def check_oneclass(self):
        data = [[0, 0], [0, 1], [1, 0], [1, 1]]
        dtype = svm.LinearOneClassData()
        model = svm.OneClassModel(dtype, nu=0.5)
        results = model.fit(data)
        for sample in data:
            print results.predict_values(sample)

        print results.predict_values([0.2, 0.2])
        print results.predict_values([2., 2.])
        print results.predict([0.2, 0.2])
        print results.predict([2., 2.])

if __name__ == '__main__':
    NumpyTest().run()
