#! /usr/bin/env python
# Last Change: Thu Nov 16 09:00 PM 2006 J

from numpy.testing import *

set_package_path()
from pyem.examples import ex1, ex2, ex3
restore_path()

# #Optional:
# set_local_path()
# # import modules that are located in the same directory as this file.
# restore_path()

class test_examples(NumpyTestCase):
    def check_ex1(self, level = 5):
        ex1()

    def check_ex2(self, level = 5):
        ex2()

    def check_ex3(self, level = 5):
        ex3()

if __name__ == "__main__":
    NumpyTest().run()
