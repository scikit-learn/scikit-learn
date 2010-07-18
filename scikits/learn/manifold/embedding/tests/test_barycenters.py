#!/usr/bin/env python

import numpy as np

from unittest import TestCase
from ..barycenters import barycenters, barycenter
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

class TestBarycenters(TestCase):
    def test_barycenters(self):
        samples = np.array((0., 0., 0.,
          1., 0., 0.,
          0., 1., 0.,
          1., 1., 0.,
          0., .5, 0.,
          .5, 0., 0.,
          1., 1., 0.5,
          )).reshape((-1,3))
        sparse = barycenters(samples, neighbors=5)
        assert(sparse.shape == (7,7))

    def test_barycenter(self):
        samples = np.array((0., 0., 0.,
          1., 0., 0.,
          0., 1., 0.,
          1., 1., 0.,
          0., .5, 0.,
          .5, 0., 0.,
          1., 1., 0.5,
          )).reshape((-1,3))
        weighs = barycenter(samples[3], samples[:3], tol=0.000001)
        assert_array_almost_equal(weighs, [-1, 1, 1], decimal=2)