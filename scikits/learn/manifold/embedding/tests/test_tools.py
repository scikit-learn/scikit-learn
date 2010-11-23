#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy
from numpy.testing import assert_array_almost_equal

from unittest import TestCase
from ..tools import create_graph

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
        g = create_graph(samples, n_neighbors=3, neigh=None,
            neigh_alternate_arguments=None)
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
        g = create_graph(samples, neigh=NewNeighbors, n_neighbors=None,
            neigh_alternate_arguments={'k' : 3})
        for l in g:
            assert(len(l) == 2)
            assert(len(l[0]) == 3)
            assert(len(l[1]) == 3)

if __name__ == "__main__":
    unittest.main()
