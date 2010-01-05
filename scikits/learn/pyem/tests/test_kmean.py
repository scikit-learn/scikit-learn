#! /usr/bin/env python
# Last Change: Thu Sep 28 01:00 PM 2006 J

import sys
from numpy.testing import *

import numpy as N

set_package_path()
from pyem.kmean import kmean
restore_path()

#Optional:
set_local_path()
# import modules that are located in the same directory as this file.
restore_path()

# Global data
X   = N.array([[3.0, 3], [4, 3], [4, 2],
        [9, 2], [5, 1], [6, 2], [9, 4], 
        [5, 2], [5, 4], [7, 4], [6, 5]])

codet1  = N.array([[3.0000, 3.0000],
        [6.2000, 4.0000], 
        [5.8000, 1.8000]])
        
codet2  = N.array([[11.0/3, 8.0/3], 
        [6.7500, 4.2500],
        [6.2500, 1.7500]])

class test_kmean(NumpyTestCase):
    def check_iter1(self, level=1):
        initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 
        code    = initc.copy()
        code1   = kmean(X, code, 1)[0]

        assert_array_almost_equal(code1, codet1)
    def check_iter2(self, level=1):
        initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 
        code    = initc.copy()
        code2   = kmean(X, code, 2)[0]

        assert_array_almost_equal(code2, codet2)

if __name__ == "__main__":
    NumpyTest().run()
