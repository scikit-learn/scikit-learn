from numpy.testing import *

# XXX remove this
import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))

import svm
import numpy as N

class test_classification(NumpyTestCase):
    def check_csvc(self):
        labels = [0, 1, 1, 2]
        samples = [[0, 0], [0, 1], [1, 0], [1, 1]]
        data = zip(labels, samples)
        dtype = svm.LinearData()
        model = svm.CSVCModel(dtype, cost=10.0)
        results = model.fit(data)
        for label, sample in data:
            assert_equal(results.predict(sample), label)
            v = results.predict_values(sample)

    def check_nusvc(self):
        labels = [0, 1, 1, 2]
        samples = [[0, 0], [0, 1], [1, 0], [1, 1]]
        data = zip(labels, samples)
        dtype = svm.LinearData()
        model = svm.NuSVCModel(dtype, nu=0.5)
        results = model.fit(data)
        for label, sample in data:
            assert_equal(results.predict(sample), label)
            v = results.predict_values(sample)

if __name__ == '__main__':
    NumpyTest().run()
