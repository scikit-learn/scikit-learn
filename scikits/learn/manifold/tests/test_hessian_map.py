#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal

from unittest import TestCase

from ...ball_tree import BallTree
from ..hessian_map import HessianMap

samples = np.array((0., 0., 0.,
  1., 0., 0.,
  0., 1., 0.,
  1., 1., 0.,
  0., .5, 0.,
  .5, 0., 0.,
  1., 1., 0.5,
  )).reshape((-1,3))

from .test_laplacian_map import close

class TestHessianMap(TestCase):
    def test_fit(self):
        np.random.seed(0)
        hessian = HessianMap(n_coords=2, n_neighbors=5)
        assert(hessian.fit(samples) == hessian)
        assert(hasattr(hessian, 'embedding_'))
        assert(hessian.embedding_.shape == (7, 2))
        neighbors_orig =\
            BallTree(samples).query(samples, k=4)[1]
        neighbors_embedding =\
            BallTree(hessian.embedding_).query(hessian.embedding_, k=4)[1]
        close(neighbors_orig, neighbors_embedding, 2)
