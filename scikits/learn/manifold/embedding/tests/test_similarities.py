#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal

from unittest import TestCase

from nose.tools import raises

from ..tools import create_neighborer
from ..similarities import LLE, HessianMap

samples = np.array((0., 0., 0.,
  1., 0., 0.,
  0., 1., 0.,
  1., 1., 0.,
  0., .5, 0.,
  .5, 0., 0.,
  1., 1., 0.5,
  )).reshape((-1,3))

from .test_similarities_mds import close

class TestLLE(TestCase):
    def test_fit(self):
        np.random.seed(0)
        lle = LLE(n_coords=2, mapping_kind=None, n_neighbors=4)
        assert(lle.fit(samples) == lle)
        assert(hasattr(lle, 'embedding_'))
        assert(lle.embedding_.shape == (7, 2))
        neighbors_orig =\
            create_neighborer(samples, n_neighbors=4).predict(samples)[1]
        neighbors_embedding =\
            create_neighborer(lle.embedding_, n_neighbors=4).predict(
                lle.embedding_)[1]
        close(neighbors_orig, neighbors_embedding, 2)

    @raises(RuntimeError)
    def test_transform_raises(self):
        np.random.seed(0)
        lle = LLE(n_coords=2, mapping_kind=None, n_neighbors=3)
        lle.fit(samples[:3])
        lle.transform(samples[0])

    def test_transform(self):
        np.random.seed(0)
        lle = LLE(n_coords=2, n_neighbors=3)
        lle.fit(samples[:3])
        mapped = lle.transform(samples)
        assert_array_almost_equal(mapped[:3], lle.embedding_, decimal=1)

class TestHessianMap(TestCase):
    def test_fit(self):
        np.random.seed(0)
        hessian = HessianMap(n_coords=2, mapping_kind=None, n_neighbors=5)
        assert(hessian.fit(samples) == hessian)
        assert(hasattr(hessian, 'embedding_'))
        assert(hessian.embedding_.shape == (7, 2))
        neighbors_orig =\
            create_neighborer(samples, n_neighbors=4).predict(samples)[1]
        neighbors_embedding =\
            create_neighborer(hessian.embedding_, n_neighbors = 4).predict(
                hessian.embedding_)[1]
        close(neighbors_orig, neighbors_embedding, 2)

    @raises(RuntimeError)
    def test_transform_raises(self):
        np.random.seed(0)
        hessian = HessianMap(n_coords=2, mapping_kind=None, n_neighbors=3)
        hessian.fit(samples[:3,:2])
        hessian.transform(samples[0,:2])

    def test_transform(self):
        np.random.seed(0)
        hessian = HessianMap(n_coords=2, n_neighbors=3)
        hessian.fit(samples[:3,:2])
        mapped = hessian.transform(samples[:,:2])
        assert_array_almost_equal(mapped[:3], hessian.embedding_, decimal=1)
