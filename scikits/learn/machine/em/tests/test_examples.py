#! /usr/bin/env python
# Last Change: Sat Jul 21 03:00 PM 2007 J

from numpy.testing import *

set_package_path()
from pyem.examples.examples import ex1, ex2, ex3, pdfestim, pdfestim1d
restore_path()

class test_examples(NumpyTestCase):
    def test_ex1(self, level = 3):
        ex1()

    def test_ex2(self, level = 3):
        ex2()

    def test_ex3(self, level = 3):
        ex3()

    def test_pdfestim(self, level = 5):
        pdfestim()

    def test_pdfestim1d(self, level = 5):
        pdfestim1d()

if __name__ == "__main__":
    NumpyTest().run()
