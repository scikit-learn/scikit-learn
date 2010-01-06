#!/usr/bin/env python

# Matthieu Brucher
# Last Change : 2007-12-10 09:52

import unittest
import numpy

from numpy.testing import *
set_package_path()
from compression import barycenters
restore_path()

class test_barycenters(unittest.TestCase):
  def test_barycenters(self):
    samples = numpy.array((0., 0., 0.,
      1., 0., 0.,
      0., 1., 0.,
      1., 1., 0.,
      0., .5, 0.,
      .5, 0., 0.,
      1., 1., 0.5,
      )).reshape((-1,3))
    sparse = barycenters(samples, neighboors=5)
    print sparse

if __name__ == "__main__":
  unittest.main()
