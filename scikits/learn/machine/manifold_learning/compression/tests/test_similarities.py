#!/usr/bin/env python

# Matthieu Brucher
# Last Change : 2007-12-10 10:14

import unittest
import numpy

from numpy.testing import *
set_package_path()
from compression import lle
restore_path()

class test_lle(unittest.TestCase):
  def test_lle(self):
    samples = numpy.array((0., 0., 0.,
      1., 0., 0.,
      0., 1., 0.,
      1., 1., 0.,
      0., .5, 0.,
      .5, 0., 0.,
      1., 1., 0.5,
      )).reshape((-1,3))
    coords = lle(samples, 2, neighboors=5)
    print coords

if __name__ == "__main__":
  unittest.main()
