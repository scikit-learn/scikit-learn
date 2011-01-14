#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal

from unittest import TestCase

from ...ball_tree import BallTree
from ..laplacian_map import LaplacianEigenmap, DiffusionMap

samples = np.array((0., 0., 0.,
  1., 0., 0.,
  0., 1., 0.,
  1., 1., 0.,
  0., .5, 0.,
  .5, 0., 0.,
  1., 1., 0.5,
  )).reshape((-1, 3))


def close(neighbor_graph_orig, neighbor_graph_estimated, value):
    for (orig, estimated) in zip(neighbor_graph_orig,
        neighbor_graph_estimated):
        assert(len(set(orig).intersection(estimated)) >= value)


class TestLaplacianEigenmap(TestCase):
    def test_fit(self):
        laplacian =\
            LaplacianEigenmap(n_coords=2, n_neighbors=4)
        assert(laplacian.fit(samples) == laplacian)
        assert(hasattr(laplacian, 'embedding_'))
        assert(laplacian.embedding_.shape == (7, 2))
        neighbors_orig =\
            BallTree(samples).query(samples, k=4)[1]
        neighbors_embedding =\
            BallTree(laplacian.embedding_).query(laplacian.embedding_, k=4)[1]
        close(neighbors_orig, neighbors_embedding, 2)


class TestDiffusionMap(TestCase):
    def test_fit(self):
        diffusion = DiffusionMap(n_coords=2, kernel_width=1)
        assert(diffusion.fit(samples) == diffusion)
        assert(hasattr(diffusion, 'embedding_'))
        assert(diffusion.embedding_.shape == (7, 2))
        neighbors_orig =\
            BallTree(samples).query(samples, k=4)[1]
        neighbors_embedding =\
            BallTree(diffusion.embedding_).query(diffusion.embedding_, k=4)[1]
        close(neighbors_orig, neighbors_embedding, 2)
