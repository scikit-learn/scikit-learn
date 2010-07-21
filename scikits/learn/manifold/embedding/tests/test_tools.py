#!/usr/bin/env python

import unittest
import numpy

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
    g = create_graph(samples, n_neighbors = 3)
    for l in g:
      assert(len(l) == 2)
      assert(len(l[0]) == 3)
      assert(len(l[1]) == 3)

if __name__ == "__main__":
  unittest.main()
  