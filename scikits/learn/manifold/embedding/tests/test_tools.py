#!/usr/bin/env python

import unittest
import numpy
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

from unittest import TestCase
from ..tools import create_graph, dist2hd

class TestCreateGraph(TestCase):
    def test_main(self):
        samples = numpy.array((0., 0., 0.,
          1., 0., 0.,
          0., 1., 0.,
          1., 1., 0.,
          0., .5, 0.,
          .5, 0., 0.,
          1., 1., 0.5,
          )).reshape((-1,3))
        g = create_graph(samples, n_neighbors = 3, neigh = None,
            neigh_alternate_arguments = None)
        for l in g:
            assert(len(l) == 2)
            assert(len(l[0]) == 3)
            assert(len(l[1]) == 3)

    def test_alternate_main(self):
        from .pseudo_neighbor import NewNeighbors
        samples = numpy.array((0., 0., 0.,
          1., 0., 0.,
          0., 1., 0.,
          1., 1., 0.,
          0., .5, 0.,
          .5, 0., 0.,
          1., 1., 0.5,
          )).reshape((-1,3))
        g = create_graph(samples, neigh = NewNeighbors, n_neighbors = None,
            neigh_alternate_arguments = {'k' : 3})
        for l in g:
            assert(len(l) == 2)
            assert(len(l[0]) == 3)
            assert(len(l[1]) == 3)

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
  