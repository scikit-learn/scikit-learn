import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

from sklearn.metrics import euclidean_distances
from sklearn.manifold import _mds as mds


def test_smacof():
    # test metric smacof using the data of "Modern Multidimensional Scaling",
    # Borg & Groenen, p 154
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    Z = np.array([[-.266, -.539],
                  [.451, .252],
                  [.016, -.238],
                  [-.200, .524]])
    X, _ = mds.smacof(sim, init=Z, n_components=2, max_iter=1, n_init=1)
    X_true = np.array([[-1.415, -2.471],
                       [1.633, 1.107],
                       [.249, -.067],
                       [-.468, 1.431]])
    assert_array_almost_equal(X, X_true, decimal=3)


def test_smacof_error():
    # Not symmetric similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim)

    # Not squared similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim)

    # init not None and not correct format:
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    Z = np.array([[-.266, -.539],
                  [.016, -.238],
                  [-.200, .524]])
    with pytest.raises(ValueError):
        mds.smacof(sim, init=Z, n_init=1)


def test_MDS():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    mds_clf = mds.MDS(metric=False, n_jobs=3, dissimilarity="precomputed")
    mds_clf.fit(sim)


def test_tau():
    sim = np.array([[0., 1., 2., 1.],
                    [1., 0., 1., 2.],
                    [2., 1., 0., 1.],
                    [1., 2., 1., 0.]])
    w = np.ones((sim.shape[0], 1))
    B = mds._tau(w, sim)
    B_true = np.array([[0.5, 0., -0.5, 0.],
                       [0., 0.5, 0., -0.5],
                       [-0.5, 0., 0.5, 0.],
                       [0., -0.5, -0., 0.5]])
    assert_array_almost_equal(B, B_true, decimal=3)

    n = 5
    sim = np.random.rand(n, n)
    w = np.ones((n, 1))
    assert_array_almost_equal(np.dot(mds._tau(w, sim), w), np.zeros((n, 1)),
                              decimal=3)


def test_tau_restriction():
    sim = np.array([[-0.05423067,  0.03647686,  0.01775381,  0.3335165],
                    [-0.01447033,  0.00349478,  0.01097555, -0.14864747],
                    [0.06870101, -0.03997165, -0.02872936, -0.18486903],
                    [-0.03937655, -0.07086947,  0.11024602,  0.22718871]])

    w = np.concatenate([np.ones(sim.shape[0] - 1), np.zeros(1)])[:, np.newaxis]
    B_proj = mds._tau(w, sim)

    w_sub = w[:3]
    sim_sub = sim[:3, :3]
    B_sub = mds._tau(w_sub, sim_sub)
    assert_array_almost_equal(B_sub, B_proj[:3, :3])


def test_MDS_transform_sample():
    data = np.array([[0., 0.],
                     [1., 0.],
                     [0., 1.]])

    n_components = 2
    clf = mds.MDS(n_components, max_iter=10000, n_init=1)

    data_embedded = clf.fit_transform(np.sqrt(data))
    data_new = np.array([[.5, .5]])
    data_embedded_new = clf.transform(data_new)

    diss_true = euclidean_distances(np.concatenate([data, data_new]))
    diss_embed = euclidean_distances(np.concatenate([data_embedded,
                                                     data_embedded_new]))

    # In-sample data constrains placement of out-of-sample data
    assert_array_almost_equal(diss_true, diss_embed, decimal=0)


def test_MDS_transform_precomputed():
    data_true = np.array([[0., 0.],
                          [0., 1.],
                          [1., 1.],
                          [1., 0.]])
    data_new = np.array([[.5, .5]])

    # Squared distance matrix for points on unit square in R^2
    sim = np.array([[0., 1., 2., 1.],
                    [1., 0., 1., 2.],
                    [2., 1., 0., 1.],
                    [1., 2., 1., 0.]])

    clf = mds.MDS(2, max_iter=100000, dissimilarity='precomputed')
    data_embed = clf.fit_transform(np.sqrt(sim))

    # Squared distances with additional center point
    sim_oos = euclidean_distances(np.concatenate([data_true, data_new]),
                                  squared=True)

    data_new_embed = clf.transform(np.sqrt(sim_oos), data_type='dissimilarity')

    data_embed_all = np.concatenate([data_embed, data_new_embed])
    sim_embed = euclidean_distances(data_embed_all, squared=True)

    assert_array_almost_equal(sim_oos[:, -1:], sim_embed[:, -1:], decimal=1)
    assert_array_almost_equal(sim_oos[-1:, :], sim_embed[-1:, :], decimal=1)
