from unittest.mock import Mock

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_almost_equal, assert_equal

from sklearn.datasets import load_digits, load_iris
from sklearn.manifold import ClassicalMDS
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
    # Testing that nonmetric MDS results in lower normalized stress compared
    # compared to metric MDS (non-regression test for issue 27028)
    X, _ = load_iris(return_X_y=True)
    sim = euclidean_distances(X)
    np.random.seed(42)
    Z = np.random.normal(size=(X.shape[0], 2))

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

    # A metric MDS solution (local minimum of the raw stress) can be rescaled to
    # decrease the stress-1 (which is returned with normalized_stress=True).
    # The optimal rescaling can be computed analytically, see Borg & Groenen,
    # Modern Multidimensional Scaling, Chapter 11.1. After rescaling, stress-1
    # becomes sqrt(s^2 / (1 + s^2)), where s is the value of stress-1 before
    # rescaling.
    stress1_rescaled = np.sqrt(stress1**2 / (1 + stress1**2))
    assert stress1_rescaled > stress2


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
        max_iter=2,
        metric_mds=False,
        init="random",
        random_state=42,
    ).fit(X)
    stress_after_2_iter = mds_est.stress_

    mds_est = mds.MDS(
        n_components=2,
        n_init=1,
        max_iter=3,
        metric_mds=False,
        init="random",
        random_state=42,
    ).fit(X)
    stress_after_3_iter = mds_est.stress_

    assert stress_after_2_iter > stress_after_3_iter


@pytest.mark.parametrize("metric_mds", [True, False])
def test_mds_recovers_true_data(metric_mds):
    X = np.array([[1, 1], [1, 4], [1, 5], [3, 3]])
    mds_est = mds.MDS(
        n_components=2,
        n_init=1,
        eps=1e-15,
        max_iter=1000,
        metric_mds=metric_mds,
        init="random",
        random_state=42,
    ).fit(X)
    stress = mds_est.stress_
    assert_allclose(stress, 0, atol=1e-6)


def test_smacof_error():
    # Not symmetric similarity matrix:
    sim = np.array([[0, 5, 9, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim, n_init=1)

    # Not squared similarity matrix:
    sim = np.array([[0, 5, 9, 4], [5, 0, 2, 2], [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim, n_init=1)

    # init not None and not correct format:
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    Z = np.array([[-0.266, -0.539], [0.016, -0.238], [-0.200, 0.524]])
    with pytest.raises(ValueError):
        mds.smacof(sim, init=Z, n_init=1)


# TODO: remove mark once loky bug is fixed:
# https://github.com/joblib/loky/issues/458
@pytest.mark.thread_unsafe
def test_MDS():
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])
    mds_clf = mds.MDS(
        metric_mds=False,
        n_jobs=3,
        n_init=3,
        metric="precomputed",
        init="random",
    )
    mds_clf.fit(sim)


# TODO(1.10): remove warning filter
@pytest.mark.filterwarnings("ignore::FutureWarning")
@pytest.mark.parametrize("k", [0.5, 1.5, 2])
def test_normed_stress(k):
    """Test that non-metric MDS normalized stress is scale-invariant."""
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    X1, stress1 = mds.smacof(sim, metric=False, max_iter=5, random_state=0)
    X2, stress2 = mds.smacof(k * sim, metric=False, max_iter=5, random_state=0)

    assert_allclose(stress1, stress2, rtol=1e-5)
    assert_allclose(X1, X2, rtol=1e-5)


# TODO(1.10): remove warning filter
@pytest.mark.filterwarnings("ignore::FutureWarning")
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


# TODO(1.10): remove warning filter
@pytest.mark.filterwarnings("ignore::FutureWarning")
@pytest.mark.parametrize("normalized_stress", [True, False])
def test_returned_stress(normalized_stress):
    # Test that the final stress corresponds to the final embedding
    # (non-regression test for issue 16846)
    X = np.array([[1, 1], [1, 4], [1, 5], [3, 3]])
    D = euclidean_distances(X)

    mds_est = mds.MDS(
        n_components=2,
        random_state=42,
        normalized_stress=normalized_stress,
    ).fit(X)

    Z = mds_est.embedding_
    stress = mds_est.stress_

    D_mds = euclidean_distances(Z)
    stress_Z = ((D_mds.ravel() - D.ravel()) ** 2).sum() / 2

    if normalized_stress:
        stress_Z = np.sqrt(stress_Z / ((D_mds.ravel() ** 2).sum() / 2))

    assert_allclose(stress, stress_Z)


# TODO(1.10): remove warning filter
@pytest.mark.filterwarnings("ignore::FutureWarning")
@pytest.mark.parametrize("metric_mds", [True, False])
def test_convergence_does_not_depend_on_scale(metric_mds):
    # Test that the number of iterations until convergence does not depend on
    # the scale of the input data
    X = np.array([[1, 1], [1, 4], [1, 5], [3, 3]])

    mds_est = mds.MDS(
        n_components=2,
        random_state=42,
        metric_mds=metric_mds,
    )

    mds_est.fit(X * 100)
    n_iter1 = mds_est.n_iter_

    mds_est.fit(X / 100)
    n_iter2 = mds_est.n_iter_

    assert_equal(n_iter1, n_iter2)


# TODO(1.9): delete this test
def test_future_warning_n_init():
    X = np.array([[1, 1], [1, 4], [1, 5], [3, 3]])
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    with pytest.warns(FutureWarning):
        mds.smacof(sim)

    with pytest.warns(FutureWarning):
        mds.MDS(init="random").fit(X)


# TODO(1.9): delete the n_init warning check
# TODO(1.10): delete this test
def test_future_warning_init_and_metric():
    X = np.array([[1, 1], [1, 4], [1, 5], [3, 3]])
    sim = np.array([[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]])

    # dissimilarity argument deprecated
    with pytest.warns(FutureWarning, match="`dissimilarity` parameter is"):
        mds.MDS(dissimilarity="precomputed", init="random", n_init=1).fit(sim)

    # metric=True deprecated
    with pytest.warns(FutureWarning, match="Use metric_mds"):
        mds.MDS(metric=True, init="random", n_init=1).fit(X)

    # metric=False deprecated
    with pytest.warns(FutureWarning, match="Use metric_mds"):
        mds.MDS(metric=False, init="random", n_init=1).fit(X)

    # default init will become classical_mds in the future
    with pytest.warns(FutureWarning, match="The default value of `init`"):
        mds.MDS(metric="euclidean", n_init=1).fit(X)

    # TODO (1.9): delete this check
    # n_init=1 will become default in the future
    with pytest.warns(FutureWarning, match="The default value of `n_init`"):
        mds.MDS(metric="euclidean", init="random").fit(X)

    # providing both metric and dissimilarity raises an error
    with pytest.raises(ValueError, match="provided both `dissimilarity`"):
        mds.MDS(
            metric="cosine", dissimilarity="euclidean", init="random", n_init=1
        ).fit(X)


# TODO(1.9): remove warning filter
@pytest.mark.filterwarnings("ignore::FutureWarning")
def test_classical_mds_init_to_mds():
    X, _ = load_iris(return_X_y=True)

    cmds = ClassicalMDS()
    Z_classical = cmds.fit_transform(X)

    mds1 = mds.MDS(init="classical_mds")
    Z1 = mds1.fit_transform(X)

    mds2 = mds.MDS(init="random")
    Z2 = mds1.fit_transform(X, init=Z_classical)

    assert_allclose(Z1, Z2)
