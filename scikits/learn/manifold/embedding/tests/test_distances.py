#!/usr/bin/env python

import unittest
import numpy

from unittest import TestCase
from ..distances import parzen, kneigh, numpy_floyd

class TestParzen(TestCase):
  def test_main(self):
    samples = numpy.array((0., 0., 0.,
      1., 0., 0.,
      0., 1., 0.,
      1., 1., 0.,
      0., .5, 0.,
      .5, 0., 0.,
      1., 1., 0.5,
      )).reshape((-1,3))
    neighbors = parzen(samples, .75)
    for i, neighlist in enumerate(neighbors):
      for neighbor in neighlist:
        d = samples[i] - samples[neighbor]
        assert( numpy.sum(d**2) < .75**2 )

class TestKNeigh(TestCase):
  def test_main(self):
    samples = numpy.array((0., 0., 0.,
      1., 0., 0.,
      0., 1., 0.,
      1., 1., 0.,
      0., .5, 0.,
      .5, 0., 0.,
      1., 1., 0.5,
      )).reshape((-1,3))
    neighbors = kneigh(samples, 2)
    for i, neighlist in enumerate(neighbors):
      assert(len(neighlist) == 2)

class TestNumpyFloyd(TestCase):
  pass

if __name__ == "__main__":
  unittest.main()
  