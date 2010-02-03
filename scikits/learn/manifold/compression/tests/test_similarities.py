#!/usr/bin/env python

# Matthieu Brucher
# Last Change : 2007-12-10 10:14

import unittest
import numpy

from unittest import TestCase
from ..similarities import LLE

class test_lle(TestCase):
  def test_lle(self):
    samples = numpy.array((0., 0., 0.,
      1., 0., 0.,
      0., 1., 0.,
      1., 1., 0.,
      0., .5, 0.,
      .5, 0., 0.,
      1., 1., 0.5,
      )).reshape((-1,3))
    coords = LLE(samples, 2, k=5)
    print coords
