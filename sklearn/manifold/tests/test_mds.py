import numpy as np
from numpy.testing import assert_array_almost_equal, assert_allclose

from sklearn.manifold import mds
from sklearn.utils.testing import assert_raises
from sklearn.utils import check_random_state


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

    assert_raises(ValueError, mds.smacof, sim)

    # Not squared similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [4, 2, 1, 0]])

    assert_raises(ValueError, mds.smacof, sim)

    # init not None and not correct format:
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    Z = np.array([[-.266, -.539],
                  [.016, -.238],
                  [-.200, .524]])
    assert_raises(ValueError, mds.smacof, sim, init=Z, n_init=1)


def test_MDS():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    mds_clf = mds.MDS(metric=False, n_jobs=3, dissimilarity="precomputed")
    mds_clf.fit(sim)


def test_normed_stress():
    # Generate a random symmetric matrix
    random_state = check_random_state(10)
    N = 100
    sim = random_state.randint(0, 2000, size=(N, N))
    sim = (sim + sim.T)/2

    # Calculate normed stress for matrix
    # and its multiplied copy
    k = 2
    stress1 = mds.smacof(sim, normalize=True, random_state=0)[1]
    stress2 = mds.smacof(k * sim, normalize=True, random_state=0)[1]

    # Normed stress should be the same for
    # values multiplied by some factor "k"
    # for a big enough matrix starting at same init
    # configuration of the embedding to initialize the
    # SMACOF algorithm..
    assert_allclose(stress1, stress2)
