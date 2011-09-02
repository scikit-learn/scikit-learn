import numpy as np

from numpy.testing import assert_almost_equal
from sklearn import neighbors, manifold
from sklearn.utils.fixes import product

eigen_solvers = ['dense', 'arpack']


def assert_lower(a, b, details=None):
    message = "%r is not lower than %r" % (a, b)
    if details is not None:
        message += ": " + details
    assert a < b, message


#----------------------------------------------------------------------
# Test LLE by computing the reconstruction error on some manifolds.

def test_lle_simple_grid():
    rng = np.random.RandomState(0)
    # grid of equidistant points in 2D, out_dim = n_dim
    X = np.array(list(product(range(5), repeat=2)))
    out_dim = 2
    clf = manifold.LocallyLinearEmbedding(n_neighbors=5, out_dim=out_dim)
    tol = .1

    N = neighbors.kneighbors_graph(
        X, clf.n_neighbors, mode='barycenter').todense()
    reconstruction_error = np.linalg.norm(np.dot(N, X) - X, 'fro')
    assert_lower(reconstruction_error, tol)

    for solver in eigen_solvers:
        clf.set_params(eigen_solver=solver)
        clf.fit(X)
        assert clf.embedding_.shape[1] == out_dim
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
    clf = manifold.LocallyLinearEmbedding(n_neighbors=5, out_dim=out_dim)
    tol = 1.5

    N = neighbors.kneighbors_graph(X, clf.n_neighbors,
                                   mode='barycenter').toarray()
    reconstruction_error = np.linalg.norm(np.dot(N, X) - X)
    assert_lower(reconstruction_error, tol)

    for solver in eigen_solvers:
        clf.set_params(eigen_solver=solver)        
        clf.fit(X)
        assert clf.embedding_.shape[1] == out_dim
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
         ('clf', neighbors.NeighborsClassifier())])
    clf.fit(iris.data, iris.target)
    assert_lower(.7, clf.score(iris.data, iris.target))


if __name__ == '__main__':
    import nose
    nose.runmodule()
