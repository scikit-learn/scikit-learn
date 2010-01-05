from numpy.testing import *

# XXX remove this
import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))

import svm
import numpy as N

class test_regression(NumpyTestCase):
    def check_epsilon_svr(self):
        y = [10., 20., 30., 40.]
        x = [[0, 0], [0, 1], [1, 0], [1, 1]]
        data = zip(y, x)
        dtype = svm.LinearData()
        model = svm.EpsilonSVRModel(dtype, cost=10.0, epsilon=0.1)
        results = model.fit(data)
        for label, sample in data:
            print results.predict(sample)

if __name__ == '__main__':
    NumpyTest().run()
