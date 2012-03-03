from itertools import product
import numpy as np
from nose.tools import assert_true

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from sklearn import neighbors, manifold
from sklearn.manifold.locally_linear import barycenter_kneighbors_graph

eigen_solvers = ['dense', 'arpack']


def assert_lower(a, b, details=None):
    message = "%r is not lower than %r" % (a, b)
    if details is not None:
        message += ": " + details
    assert a < b, message


#----------------------------------------------------------------------
# Test utility routines
def test_barycenter_kneighbors_graph():
    X = np.array([[0, 1], [1.01, 1.], [2, 0]])

    A = barycenter_kneighbors_graph(X, 1)
    assert_array_almost_equal(
        A.todense(),
        [[0.,  1.,  0.],
         [1.,  0.,  0.],
         [0.,  1.,  0.]])

    A = barycenter_kneighbors_graph(X, 2)
    # check that columns sum to one
    assert_array_almost_equal(np.sum(A.todense(), 1), np.ones((3, 1)))
    pred = np.dot(A.todense(), X)
    assert_true(np.linalg.norm(pred - X) / X.shape[0] < 1)


#----------------------------------------------------------------------
# Test LLE by computing the reconstruction error on some manifolds.

def test_lle_simple_grid():
    rng = np.random.RandomState(0)
    # grid of equidistant points in 2D, out_dim = n_dim
    X = np.array(list(product(range(5), repeat=2)))
    out_dim = 2
    clf = manifold.LocallyLinearEmbedding(n_neighbors=5, out_dim=out_dim)
    tol = .1

    N = barycenter_kneighbors_graph(X, clf.n_neighbors).todense()
    reconstruction_error = np.linalg.norm(np.dot(N, X) - X, 'fro')
    assert_lower(reconstruction_error, tol)

    for solver in eigen_solvers:
        clf.set_params(eigen_solver=solver)
        clf.fit(X)
        assert_true(clf.embedding_.shape[1] == out_dim)
        reconstruction_error = np.linalg.norm(
            np.dot(N, clf.embedding_) - clf.embedding_, 'fro') ** 2
        # FIXME: ARPACK fails this test ...
        if solver != 'arpack':
            assert_lower(reconstruction_error, tol)
            assert_almost_equal(clf.reconstruction_error_,
                                reconstruction_error, decimal=4)

    # re-embed a noisy version of X using the transform method
    noise = rng.randn(*X.shape) / 100
    X_reembedded = clf.transform(X + noise)
    assert_lower(np.linalg.norm(X_reembedded - clf.embedding_), tol)


def test_lle_manifold():
    # similar test on a slightly more complex manifold
    X = np.array(list(product(range(20), repeat=2)))
    X = np.c_[X, X[:, 0] ** 2 / 20]
    out_dim = 2
    clf = manifold.LocallyLinearEmbedding(n_neighbors=5, out_dim=out_dim,
                                          random_state=0)
    tol = 1.5

    N = barycenter_kneighbors_graph(X, clf.n_neighbors).toarray()
    reconstruction_error = np.linalg.norm(np.dot(N, X) - X)
    assert_lower(reconstruction_error, tol)

    for solver in eigen_solvers:
        clf.set_params(eigen_solver=solver)
        clf.fit(X)
        assert_true(clf.embedding_.shape[1] == out_dim)
        reconstruction_error = np.linalg.norm(
            np.dot(N, clf.embedding_) - clf.embedding_, 'fro') ** 2
        details = "solver: " + solver
        assert_lower(reconstruction_error, tol, details=details)
        assert_lower(np.abs(clf.reconstruction_error_ - reconstruction_error),
                     tol * reconstruction_error, details=details)


def test_pipeline():
    # check that LocallyLinearEmbedding works fine as a Pipeline
    from sklearn import pipeline, datasets
    iris = datasets.load_iris()
    clf = pipeline.Pipeline(
        [('filter', manifold.LocallyLinearEmbedding()),
         ('clf', neighbors.KNeighborsClassifier())])
    clf.fit(iris.data, iris.target)
    assert_lower(.7, clf.score(iris.data, iris.target))


# Test the error raised when the weight matrix is singular
def test_singular_matrix():
    from nose.tools import assert_raises
    M = np.ones((10, 3))

    assert_raises(ValueError, manifold.locally_linear_embedding,
                  M, 2, 1, method='standard', eigen_solver='arpack')


if __name__ == '__main__':
    import nose
    nose.runmodule()
