#!/usr/bin/env python

import unittest
import numpy
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

from unittest import TestCase
from nose.tools import raises

from ..builder import builder
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

class TestBuilder(TestCase):
    def test_no_mapping(self):
        test = builder(None, kind = None)
        assert(test is None)

    def test_default_construction(self):
        embedding = Embedding()
        embedding.X_ = samples
        embedding.embedding_ = samples
        test = builder(embedding, n_neighbors=3)
        assert(isinstance(test, Barycenter))

    def test_construction_class(self):
        embedding = Embedding()
        embedding.X_ = samples
        embedding.embedding_ = samples
        test = builder(embedding, kind = Barycenter, n_neighbors=3)
        assert(isinstance(test, Barycenter))

    def test_construction_instance(self):
        test = Barycenter(n_neighbors=3, tol=1e-6)
        embedding = Embedding()
        embedding.X_ = samples
        embedding.embedding_ = samples
        test = builder(embedding, kind = test, n_neighbors=3)
        assert(isinstance(test, Barycenter))
    
    @raises(RuntimeError)
    def test_construction_exception(self):
        embedding = Embedding()
        embedding.X_ = samples
        embedding.embedding_ = samples
        test = builder(embedding, kind = embedding, n_neighbors=3)
        assert(isinstance(test, Barycenter))

if __name__ == "__main__":
  unittest.main()
  