#! /usr/bin/env python
# Last Change: Sun Sep 07 04:00 PM 2008 J

from unittest import TestCase

from scikits.learn.machine.em.examples.examples import ex1, ex2, ex3, pdfestim, pdfestim1d

class test_examples(TestCase):
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
