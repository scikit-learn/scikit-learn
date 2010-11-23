#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal

from unittest import TestCase

from nose.tools import raises

from ..tools import create_neighborer
from ..similarities_mds import LaplacianEigenmap, DiffusionMap

samples = np.array((0., 0., 0.,
  1., 0., 0.,
  0., 1., 0.,
  1., 1., 0.,
  0., .5, 0.,
  .5, 0., 0.,
  1., 1., 0.5,
  )).reshape((-1,3))

def close(neighbor_graph_orig, neighbor_graph_estimated, value):
    for (orig, estimated) in zip(neighbor_graph_orig, neighbor_graph_estimated):
        assert(len(set(orig).intersection(estimated)) >= value)

class TestLaplacianEigenmap(TestCase):
    def test_fit(self):
        laplacian =\
            LaplacianEigenmap(n_coords=2, mapping_kind=None, n_neighbors=4)
        assert(laplacian.fit(samples) == laplacian)
        assert(hasattr(laplacian, 'embedding_'))
        assert(laplacian.embedding_.shape == (7, 2))
        neighbors_orig =\
            create_neighborer(samples, n_neighbors=4).predict(samples)[1]
        neighbors_embedding =\
            create_neighborer(laplacian.embedding_, n_neighbors=4).predict(
                laplacian.embedding_)[1]
        close(neighbors_orig, neighbors_embedding, 2)

    @raises(RuntimeError)
    def test_transform_raises(self):
        laplacian =\
             LaplacianEigenmap(n_coords=2, mapping_kind=None, n_neighbors=4)
        laplacian.fit(samples[:4, :2])
        laplacian.transform(samples[0])

    def test_transform(self):
        laplacian = LaplacianEigenmap(n_coords=2, n_neighbors=4)
        laplacian.fit(samples[:4, :2])
        mapped = laplacian.transform(samples[:, :2])
        assert_array_almost_equal(mapped[:4, :2],
            laplacian.embedding_, decimal=1)

class TestDiffusionMap(TestCase):
    def test_fit(self):
        diffusion = DiffusionMap(n_coords=2, mapping_kind=None, kernel_width=1)
        assert(diffusion.fit(samples) == diffusion)
        assert(hasattr(diffusion, 'embedding_'))
        assert(diffusion.embedding_.shape == (7, 2))
        neighbors_orig =\
            create_neighborer(samples, n_neighbors=4).predict(samples)[1]
        neighbors_embedding =\
            create_neighborer(diffusion.embedding_, n_neighbors=4).predict(
                diffusion.embedding_)[1]
        close(neighbors_orig, neighbors_embedding, 2)

    @raises(RuntimeError)
    def test_transform_raises(self):
        diffusion = DiffusionMap(n_coords=2, mapping_kind=None, n_neighbors=3)
        diffusion.fit(samples[:3, :2])
        diffusion.transform(samples[0, :2])

    def test_transform(self):
        diffusion = DiffusionMap(n_coords=2, n_neighbors=3)
        diffusion.fit(samples[:3, :2])
        mapped = diffusion.transform(samples[:, :2])
        assert_array_almost_equal(mapped[:3, :2],
            diffusion.embedding_, decimal=1)
