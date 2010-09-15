#!/usr/bin/env python

import unittest
import numpy
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

from unittest import TestCase
from ..barycenter import Barycenter

samples = numpy.array((0., 0., 0.,
  1., 0., 0.,
  0., 1., 0.,
  1., 1., 0.,
  0., .5, 0.,
  .5, 0., 0.,
  1., 1., 0.5,
  )).reshape((-1,3))

test_data = ((.4, .5, 0), (.5, .4, 0), (.5, .5, 0), )

class Embedding(object):
    pass

class TestBarycenter(TestCase):
    def test_create(self):
        test = Barycenter(n_neighbors=3)

    def test_fit(self):
        test = Barycenter(n_neighbors=3)
        embedding = Embedding()
        embedding.X_ = samples
        embedding.embedding_ = samples
        r = test.fit(embedding)
        assert(r is test)

    def test_transform_1D(self):
        test = Barycenter(n_neighbors=3, tol=1e-6)
        embedding = Embedding()
        embedding.X_ = samples
        embedding.embedding_ = samples
        test.fit(embedding)
        embedded = test.transform(test_data[2])
        assert_array_almost_equal(embedded, numpy.atleast_2d(test_data[2]),
            decimal=5)

    def test_transform_2D(self):
        test = Barycenter(n_neighbors=3, tol=1e-6)
        embedding = Embedding()
        embedding.X_ = samples
        embedding.embedding_ = samples
        test.fit(embedding)
        embedded = test.transform(test_data)
        assert_array_almost_equal(embedded, test_data, decimal=5)

if __name__ == "__main__":
  unittest.main()
  