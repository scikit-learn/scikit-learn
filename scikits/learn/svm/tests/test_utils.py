from numpy.testing import *

# XXX remove this
import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))

import svm.utils as utils
import numpy as N

class test_util(NumpyTestCase):
    def check_addressof_array(self):
        a = N.array([1.])
        addr = utils.addressof_array(a)
        self.assert_(addr != 0)

if __name__ == '__main__':
    NumpyTest().run()
