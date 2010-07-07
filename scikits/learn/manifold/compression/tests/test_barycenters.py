#!/usr/bin/env python

# Matthieu Brucher
# Last Change : 2007-12-10 09:52

import numpy as np

from unittest import TestCase
from ..barycenters import barycenters

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
    sparse = barycenters(samples, k=5)
    print sparse
