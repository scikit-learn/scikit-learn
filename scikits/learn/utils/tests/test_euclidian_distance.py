#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy
from numpy.testing import assert_array_almost_equal

from unittest import TestCase
from ..euclidian_distance import dist2hd

class TestDist2HD(TestCase):
    def test_dist2hd(self):
        samples = numpy.array((0., 0., 0.,
          1., 0., 0.,
          0., 1., 0.,
          1., 1., 0.,
          0., .5, 0.,
          .5, 0., 0.,
          1., 1., 0.5,
          )).reshape((-1,3))
        distances = dist2hd(samples, numpy.array(((0., 0., 0.), )))
        assert_array_almost_equal(distances**2,
            numpy.array(((0, 1, 1, 2, .25, .25, 2.25), )).T)

if __name__ == "__main__":
    unittest.main()
