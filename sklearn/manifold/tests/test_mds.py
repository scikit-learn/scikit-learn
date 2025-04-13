from unittest.mock import Mock

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_almost_equal

from sklearn.datasets import load_digits
from sklearn.manifold import _mds as mds
from sklearn.metrics import euclidean_distances


def test_smacof():
    # test metric smacof using the data of "Modern Multidimensional Scaling",
    # Borg & Groenen, p 154
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])
    Z = np.array([[-0.266, -0.539], [0.451, 0.252], [0.016, -0.238], [-0.200, 0.524]])
    X, _ = mds.smacof(sim, init=Z, n_components=2, max_iter=1, n_init=1)
    X_true = np.array(
        [[-1.415, -2.471], [1.633, 1.107], [0.249, -0.067], [-0.468, 1.431]]
    )
    assert_array_almost_equal(X, X_true, decimal=3)


def test_nonmetric_lower_normalized_stress():
    # Testing that nonmetric MDS results in lower normalized stess compared
    # compared to metric MDS (non-regression test for issue 27028)
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])
    Z = np.array([[-0.266, -0.539], [0.451, 0.252], [0.016, -0.238], [-0.200, 0.524]])

    _, stress1 = mds.smacof(
        sim, init=Z, n_components=2, max_iter=1000, n_init=1, normalized_stress=True
    )

    _, stress2 = mds.smacof(
        sim,
        init=Z,
        n_components=2,
        max_iter=1000,
        n_init=1,
        normalized_stress=True,
        metric=False,
    )
    assert stress1 > stress2


def test_nonmetric_mds_optimization():
    # Test that stress is decreasing during nonmetric MDS optimization
    # (non-regression test for issue 27028)
    X, _ = load_digits(return_X_y=True)
    rng = np.random.default_rng(seed=42)
    ind_subset = rng.choice(len(X), size=200, replace=False)
    X = X[ind_subset]

    mds_est = mds.MDS(
        n_components=2,
        n_init=1,
        eps=1e-15,
        max_iter=2,
        metric=False,
        random_state=42,
    ).fit(X)
    stress_after_2_iter = mds_est.stress_

    mds_est = mds.MDS(
        n_components=2,
        n_init=1,
        eps=1e-15,
        max_iter=3,
        metric=False,
        random_state=42,
    ).fit(X)
    stress_after_3_iter = mds_est.stress_

    assert stress_after_2_iter > stress_after_3_iter


@pytest.mark.parametrize("metric", [True, False])
def test_mds_recovers_true_data(metric):
    X = np.array([[1, 1], [1, 4], [1, 5], [3, 3]])
    mds_est = mds.MDS(
        n_components=2,
        n_init=1,
        eps=1e-15,
        max_iter=1000,
        metric=metric,
        random_state=42,
    ).fit(X)
    stress = mds_est.stress_
    assert_allclose(stress, 0, atol=1e-10)


def test_smacof_error():
    # Not symmetric similarity matrix:
    sim = np.array([[0, 5, 9, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim)

    # Not squared similarity matrix:
    sim = np.array([[0, 5, 9, 4], [5, 0, 2, 2], [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim)

    # init not None and not correct format:
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    Z = np.array([[-0.266, -0.539], [0.016, -0.238], [-0.200, 0.524]])
    with pytest.raises(ValueError):
        mds.smacof(sim, init=Z, n_init=1)


def test_MDS():
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])
    mds_clf = mds.MDS(metric=False, n_jobs=3, dissimilarity="precomputed")
    mds_clf.fit(sim)


@pytest.mark.parametrize("k", [0.5, 1.5, 2])
def test_normed_stress(k):
    """Test that non-metric MDS normalized stress is scale-invariant."""
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    X1, stress1 = mds.smacof(sim, metric=False, max_iter=5, random_state=0)
    X2, stress2 = mds.smacof(k * sim, metric=False, max_iter=5, random_state=0)

    assert_allclose(stress1, stress2, rtol=1e-5)
    assert_allclose(X1, X2, rtol=1e-5)


@pytest.mark.parametrize("metric", [True, False])
def test_normalized_stress_auto(metric, monkeypatch):
    rng = np.random.RandomState(0)
    X = rng.randn(4, 3)
    dist = euclidean_distances(X)

    mock = Mock(side_effect=mds._smacof_single)
    monkeypatch.setattr("sklearn.manifold._mds._smacof_single", mock)

    est = mds.MDS(metric=metric, normalized_stress="auto", random_state=rng)
    est.fit_transform(X)
    assert mock.call_args[1]["normalized_stress"] != metric

    mds.smacof(dist, metric=metric, normalized_stress="auto", random_state=rng)
    assert mock.call_args[1]["normalized_stress"] != metric


def test_isotonic_outofbounds():
    # This particular configuration can trigger out of bounds error
    # in the isotonic regression (non-regression test for issue 26999)
    dis = np.array(
        [
            [0.0, 1.732050807568877, 1.7320508075688772],
            [1.732050807568877, 0.0, 6.661338147750939e-16],
            [1.7320508075688772, 6.661338147750939e-16, 0.0],
        ]
    )
    init = np.array(
        [
            [0.08665881585055124, 0.7939114643387546],
            [0.9959834154297658, 0.7555546025640025],
            [0.8766008278401566, 0.4227358815811242],
        ]
    )
    mds.smacof(dis, init=init, metric=False, n_init=1)


def test_returned_stress():
    # Test that the final stress corresponds to the final embedding
    # (non-regression test for issue 16846)
    X = np.array([[1, 1], [1, 4], [1, 5], [3, 3]])
    D = euclidean_distances(X)

    mds_est = mds.MDS(n_components=2, random_state=42).fit(X)
    Z = mds_est.embedding_
    stress = mds_est.stress_

    D_mds = euclidean_distances(Z)
    stress_Z = ((D_mds.ravel() - D.ravel()) ** 2).sum() / 2

    assert_allclose(stress, stress_Z)
