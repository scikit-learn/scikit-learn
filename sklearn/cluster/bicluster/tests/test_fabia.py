# -*- coding: utf-8 -*-
"""
Tests for the FABIA biclustering algorithm.
"""

import numpy as np

from scipy.sparse import csr_matrix

from sklearn.cluster.bicluster import FabiaBiclustering

from sklearn.metrics import consensus_score

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises

from sklearn.datasets import make_fabia_biclusters

from nose import SkipTest


def test_onecluster_nonoise():
    k, n, m = 1, 100, 100
    (Xn, X, zc, lc) = make_fabia_biclusters((n, m), k, (5, 5), (5, 5), 3,
                                            1e-9, 2.0, 1.0, 1e-9, 3.0, 1.0)
    model = FabiaBiclustering(1, 200, random_state=31)
    model.fit(X)
    s = consensus_score(model.biclusters_, (zc, lc))
    # score is 1 almost all the time,
    # but it does happens that we miss one entry
    assert_equal(s >= 0.8, True)


def test_onecluster_lightnoise():
    raise SkipTest("Skip test because of its long runtime.")
    k, n, m = 1, 100, 100
    (Xn, X, zc, lc) = make_fabia_biclusters((n, m), k, (5, 5), (5, 5), 1e-2,
                                            1e-9, 2.0, 1.0, 1e-9, 3.0, 1.0)
    model = FabiaBiclustering(1, 300, random_state=0)
    model.fit(Xn)
    s = consensus_score(model.biclusters_, (zc, lc))
    # score is 1 almost all the time,
    # but it does happens that we miss one entry
    assert_equal(s >= 0.8, True)


def test_paper():
    ''' Test FABIA on simulated datasets from the original publication.'''
    raise SkipTest("Skip FABIA test_paper() because of its long runtime.")
    k, n, m = 10, 100, 1000
    (Xn, X, zc, lc) = make_fabia_biclusters((n, m), k, (5, 25), (10, 210), 3.0,
                                            0.2, 2.0, 1.0, 0.2, 3.0, 1.0)
    model = FabiaBiclustering(n_clusters=13, alpha=0.01, scale=False)
    model.fit(Xn)
    s = consensus_score(model.biclusters_, (zc, lc))
    # the current R implementation manages an average score of ~0.5,
    # so expecting s > 0.4 is reasonable
    assert_equal(s >= 0.4, True)


def test_output_scaling():
    ''' After a FABIA run, the rows of Z should have variance of roughly 1.'''
    k, n, m = 10, 100, 1000
    (Xn, X, zc, lc) = make_fabia_biclusters((n, m), k, (5, 25), (10, 210), 3.0,
                                            0.2, 2.0, 1.0, 0.2, 3.0, 1.0)
    model = FabiaBiclustering(n_clusters=13, alpha=0.01, n_iter=2, scale=False)
    model.fit(Xn)
    vz = model.Z_.var(0)
    assert(np.allclose(vz, np.ones(vz.shape[0]), 0.25))  # very roughly


def test_parameter_validation():
    data = np.arange(25).reshape((5, 5))

    model = FabiaBiclustering(n_clusters=-1)
    assert_raises(ValueError, model.fit, data)

    model = FabiaBiclustering(n_iter=-1)
    assert_raises(ValueError, model.fit, data)

    model = FabiaBiclustering(alpha=1.05)
    assert_raises(ValueError, model.fit, data)

    model = FabiaBiclustering(spz=-1)
    assert_raises(ValueError, model.fit, data)

    model = FabiaBiclustering(spl=-1)
    assert_raises(ValueError, model.fit, data)

    model = FabiaBiclustering(eps=-1)
    assert_raises(ValueError, model.fit, data)

    model = FabiaBiclustering(rescale_l=100)
    assert_raises(ValueError, model.fit, data)

    model = FabiaBiclustering(scale="no")
    assert_raises(ValueError, model.fit, data)

    # no scaling allowed when input is sparse
    model = FabiaBiclustering(scale=True)
    Xs = csr_matrix([[0, 1, 1], [1, 1, 0], [0, 0, 0]], dtype=np.double)
    assert_raises(ValueError, model.fit, Xs)


def test_sparse():
    # Just tests that nothing throws, we'd have to run much
    # longer to obtain meaningful results.
    X = csr_matrix([[0, 1, 1], [1, 1, 0], [0, 0, 0]], dtype=np.double)
    model = FabiaBiclustering(n_clusters=1, n_iter=1, scale=False)
    model.fit(X)


def test_rescale_l():
    # Just tests that nothing throws, we'd have to run much
    # longer to obtain meaningful results.
    k, n, m = 10, 100, 1000
    (Xn, X, zc, lc) = make_fabia_biclusters((n, m), k, (5, 25), (10, 210), 3.0,
                                            0.2, 2.0, 1.0, 0.2, 3.0, 1.0)
    model = FabiaBiclustering(n_clusters=1, n_iter=1, rescale_l=True)
    model.fit(Xn)
