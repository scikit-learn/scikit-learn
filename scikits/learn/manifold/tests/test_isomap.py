import numpy as np

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from scikits.learn import neighbors, manifold, datasets, preprocessing
from scikits.learn.utils.fixes import product

eigen_solvers = ['auto', 'dense', 'arpack']
path_methods = ['auto', 'FW', 'D']


def assert_lower(a, b, details=None):
    message = "%r is not lower than %r" % (a, b)
    if details is not None:
        message += ": " + details
    assert a < b, message


def test_isomap_simple_grid():
    # Isomap should preserve distances when all neighbors are used
    N_per_side = 5
    Npts = N_per_side ** 2
    n_neighbors = Npts - 1

    # grid of equidistant points in 2D, out_dim = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)))

    # distances from each point to all others
    G = neighbors.kneighbors_graph(X, n_neighbors,
                                   mode='distance').toarray()

    for eigen_solver in eigen_solvers:
        for path_method in path_methods:
            clf = manifold.Isomap(n_neighbors=n_neighbors, out_dim=2,
                                  eigen_solver=eigen_solver,
                                  path_method=path_method)
            clf.fit(X)

            G_iso = neighbors.kneighbors_graph(clf.embedding_,
                                               n_neighbors,
                                               mode='distance').toarray()
            assert_array_almost_equal(G, G_iso)


def test_isomap_reconstruction_error():
    # Same setup as in test_isomap_simple_grid, with an added dimension
    N_per_side = 5
    Npts = N_per_side ** 2
    n_neighbors = Npts - 1

    # grid of equidistant points in 2D, out_dim = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)))

    # add noise in a third dimension
    rng = np.random.RandomState(0)
    noise = 0.1 * rng.randn(Npts, 1)
    X = np.concatenate((X, noise), 1)

    # compute input kernel
    G = neighbors.kneighbors_graph(X, n_neighbors,
                                   mode='distance').toarray()

    centerer = preprocessing.KernelCenterer()
    K = centerer.fit_transform(-0.5 * G ** 2)

    for eigen_solver in eigen_solvers:
        for path_method in path_methods:
            clf = manifold.Isomap(n_neighbors=n_neighbors, out_dim=2,
                                  eigen_solver=eigen_solver,
                                  path_method=path_method)
            clf.fit(X)

            # compute output kernel
            G_iso = neighbors.kneighbors_graph(clf.embedding_,
                                               n_neighbors,
                                               mode='distance').toarray()

            K_iso = centerer.fit_transform(-0.5 * G_iso ** 2)

            # make sure error agrees
            reconstruction_error = np.linalg.norm(K - K_iso) / Npts
            assert_almost_equal(reconstruction_error,
                                clf.reconstruction_error())


def test_transform():
    N = 200
    k = 10
    eps = 0.01

    # Create S-curve dataset
    X, y = datasets.samples_generator.s_curve(N)

    # Compute isomap embedding
    iso = manifold.Isomap(k, 2)
    X_iso = iso.fit_transform(X)

    # Re-embed a noisy version of the points
    rng = np.random.RandomState(0)
    noise = eps * rng.randn(*X.shape)
    X_iso2 = iso.transform(X + noise)

    # Make sure the rms error on re-embedding is comparable to eps
    assert np.sqrt(np.mean((X_iso - X_iso2) ** 2)) < 2 * eps


def test_pipeline():
    # check that LocallyLinearEmbedding works fine as a Pipeline
    from scikits.learn import pipeline, datasets
    iris = datasets.load_iris()
    clf = pipeline.Pipeline(
        [('filter', manifold.Isomap()),
         ('clf', neighbors.NeighborsClassifier())])
    clf.fit(iris.data, iris.target)
    assert_lower(.7, clf.score(iris.data, iris.target))


if __name__ == '__main__':
    import nose
    nose.runmodule()
