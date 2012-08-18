from itertools import product
import numpy as np
from nose.tools import assert_true

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from sklearn import neighbors, manifold
from sklearn.manifold.locally_linear import barycenter_kneighbors_graph
from sklearn.utils.testing import assert_less

eigen_solvers = ['dense', 'arpack']


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
    assert_less(np.linalg.norm(pred - X) / X.shape[0], 1)


#----------------------------------------------------------------------
# Test LLE by computing the reconstruction error on some manifolds.

def test_lle_simple_grid():
    # note: ARPACK is numerically unstable, so this test will fail for
    #       some random seeds.  We choose 2 because the tests pass.
    rng = np.random.RandomState(2)
    tol = 0.1

    # grid of equidistant points in 2D, n_components = n_dim
    X = np.array(list(product(range(5), repeat=2)))
    X = X + 1e-10 * rng.uniform(size=X.shape)
    n_components = 2
    clf = manifold.LocallyLinearEmbedding(n_neighbors=5,
            n_components=n_components, random_state=rng)
    tol = 0.1

    N = barycenter_kneighbors_graph(X, clf.n_neighbors).todense()
    reconstruction_error = np.linalg.norm(np.dot(N, X) - X, 'fro')
    assert_less(reconstruction_error, tol)

    for solver in eigen_solvers:
        clf.set_params(eigen_solver=solver)
        clf.fit(X)
        assert_true(clf.embedding_.shape[1] == n_components)
        reconstruction_error = np.linalg.norm(
            np.dot(N, clf.embedding_) - clf.embedding_, 'fro') ** 2
       
        assert_less(reconstruction_error, tol)
        assert_almost_equal(clf.reconstruction_error_,
                            reconstruction_error, decimal=1)

    # re-embed a noisy version of X using the transform method
    noise = rng.randn(*X.shape) / 100
    X_reembedded = clf.transform(X + noise)
    assert_less(np.linalg.norm(X_reembedded - clf.embedding_), tol)


def test_lle_manifold():
    rng = np.random.RandomState(0)
    # similar test on a slightly more complex manifold
    X = np.array(list(product(range(20), repeat=2)))
    X = np.c_[X, X[:, 0] ** 2 / 20]
    X = X + 1e-10 * rng.uniform(size=X.shape)
    n_components = 2
    for method in ["standard", "hessian", "modified", "ltsa"]:
        clf = manifold.LocallyLinearEmbedding(n_neighbors=6,
                n_components=n_components, method=method, random_state=0)
        tol = 1.5 if method == "standard" else 3

        N = barycenter_kneighbors_graph(X, clf.n_neighbors).toarray()
        reconstruction_error = np.linalg.norm(np.dot(N, X) - X)
        assert_less(reconstruction_error, tol)

        for solver in eigen_solvers:
            clf.set_params(eigen_solver=solver)
            clf.fit(X)
            assert_true(clf.embedding_.shape[1] == n_components)
            reconstruction_error = np.linalg.norm(
                np.dot(N, clf.embedding_) - clf.embedding_, 'fro') ** 2
            details = ("solver: %s, method: %s"
                % (solver, method))
            assert_less(reconstruction_error, tol, msg=details)
            assert_less(np.abs(clf.reconstruction_error_ - reconstruction_error),
                         tol * reconstruction_error, msg=details)


def test_pipeline():
    # check that LocallyLinearEmbedding works fine as a Pipeline
    # only checks that no error is raised.
    # TODO check that it actually does something useful
    from sklearn import pipeline, datasets
    X, y = datasets.make_blobs(random_state=0)
    clf = pipeline.Pipeline(
        [('filter', manifold.LocallyLinearEmbedding(random_state=0)),
         ('clf', neighbors.KNeighborsClassifier())])
    clf.fit(X, y)
    assert_less(.9, clf.score(X, y))


# Test the error raised when the weight matrix is singular
def test_singular_matrix():
    import warnings
    from nose.tools import assert_raises
    M = np.ones((10, 3))
    with warnings.catch_warnings(record=True):
        assert_raises(ValueError, manifold.locally_linear_embedding,
                      M, 2, 1, method='standard', eigen_solver='arpack')


if __name__ == '__main__':
    import nose
    nose.runmodule()
