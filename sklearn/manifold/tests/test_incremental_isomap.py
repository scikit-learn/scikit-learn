from itertools import product
import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from sklearn import datasets
from sklearn import manifold
from sklearn import neighbors
from sklearn import pipeline
from sklearn import preprocessing
from sklearn.utils.testing import assert_less

eigen_solvers = ['auto', 'dense', 'arpack']
path_methods = ['auto', 'FW', 'D']


def test_isomap_simple_grid():
    # Isomap should preserve distances when all neighbors are used
    N_per_side = 5
    Npts = N_per_side ** 2
    n_neighbors = Npts - 1

    # grid of equidistant points in 2D, n_components = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)))

    # distances from each point to all others
    G = neighbors.kneighbors_graph(X, n_neighbors,
                                   mode='distance').toarray()

    for eigen_solver in eigen_solvers:
        for path_method in path_methods:
            clf = manifold.IncrementalIsomap(n_neighbors=n_neighbors, n_components=2,
                                  eigen_solver=eigen_solver,
                                  path_method=path_method)
            clf.fit(X)
            embedding_unnormed = (clf.embedding_[:clf.n, :] / np.sqrt(clf.n) *
                                    np.sqrt(clf.kernel_pca_.lambdas_))  # convert scale back
            G_iso = neighbors.kneighbors_graph(embedding_unnormed,
                                               n_neighbors,
                                               mode='distance').toarray()
            assert_array_almost_equal(G, G_iso)


def test_isomap_incremental_grid():
    # Isomap on full dataset should be equal to incremental isomap
    # after full dataset has been processed
    N_per_side = 10
    Npts = N_per_side ** 2
    n_neighbors = 10

    # grid of equidistant points in 2D, n_components = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)), np.float)
    X += np.random.randn(X.shape[0], X.shape[1])

    # distances from each point to all others
    clf_iso = manifold.Isomap(n_neighbors=n_neighbors, n_components=2)
    iso_full = clf_iso.fit_transform(X)

    n_train = 20
    i_train = np.arange(n_train)  # select first n_train samples for training
    i_train.sort()
    i_test = np.setdiff1d(range(Npts), i_train)
    X_train = X[i_train]
    X_test = X[i_test]
    clf = manifold.IncrementalIsomap(n_neighbors=n_neighbors, n_components=2,
                                         n_max_samples=Npts)
    clf.fit(X_train)

    for i in range(len(i_test)):
        clf.partial_fit_transform(X_test[i, :].reshape(1, -1))

    iso_incremental = (clf.transform(X) / np.sqrt(clf.n) *
                                    np.sqrt(clf.kernel_pca_.lambdas_))
    #assert_array_almost_equal(clf_iso.kng_.todense(), clf.kng_.todense())  # only possible if both clf are of class IncrementalIsomap
    assert_array_almost_equal(iso_full, np.sign(iso_full[0, :]) *  # assure same signs for both embeddings
                                        np.sign(iso_incremental[0, :]) *
                                        iso_incremental)


def test_isomap_incremental_extend():
    # Same as above, only n_max_samples is set to smaller values and must be extended automatically
    N_per_side = 10
    Npts = N_per_side ** 2
    n_neighbors = 10

    # grid of equidistant points in 2D, n_components = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)), np.float)
    X += np.random.randn(X.shape[0], X.shape[1])

    # distances from each point to all others
    clf_iso = manifold.Isomap(n_neighbors=n_neighbors, n_components=2)
    iso_full = clf_iso.fit_transform(X)

    n_train = 20
    i_train = np.arange(n_train)  # select first n_train samples for training
    i_train.sort()
    i_test = np.setdiff1d(range(Npts), i_train)
    X_train = X[i_train]
    X_test = X[i_test]
    clf = manifold.IncrementalIsomap(n_neighbors=n_neighbors, n_components=2,
                                     n_max_samples=n_train, overflow_mode="extend")
    clf.fit(X_train)

    for i in range(len(i_test)):
        clf.partial_fit(X_test[i, :].reshape(1, -1))
    iso_incremental = (clf.transform(X) / np.sqrt(clf.n) *
                                    np.sqrt(clf.kernel_pca_.lambdas_))
    assert_array_almost_equal(iso_full, np.sign(iso_full[0, :]) *  # assure same signs for both embeddings
                                        np.sign(iso_incremental[0, :]) *
                                        iso_incremental)


def test_isomap_incremental_embed():
    # Same as above, only n_max_samples is set to smaller values remaining 
    # values are just embedded with oose
    N_per_side = 10
    Npts = N_per_side ** 2
    n_neighbors = 10

    # grid of equidistant points in 2D, n_components = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)), np.float)
    X += np.random.randn(X.shape[0], X.shape[1])

    # select training data
    n_train = 20
    i_train = np.arange(n_train)  # select first n_train samples for training
    i_train.sort()
    i_test = np.setdiff1d(range(Npts), i_train)
    X_train = X[i_train]
    X_test = X[i_test]

    # Full isomap
    clf_iso = manifold.Isomap(n_neighbors=n_neighbors, n_components=2)
    clf_iso.fit(X_train)
    iso_full = clf_iso.transform(X_test)

    clf = manifold.IncrementalIsomap(n_neighbors=n_neighbors, n_components=2,
                                     n_max_samples=n_train,
                                     overflow_mode="embed")
    clf.fit(X_train)
    iso_incremental = np.zeros((len(i_test), clf.n_components), dtype=np.float)
    for i in range(len(i_test)):
        iso_incremental[i, :] = clf.partial_fit_transform(X_test[i, :].reshape(1, -1))

    iso_incremental *= 1.0 / np.sqrt(clf.n) * np.sqrt(clf.kernel_pca_.lambdas_)

    assert_array_almost_equal(iso_full, np.sign(iso_full[0, :]) *  # assure same signs for both embeddings
                                        np.sign(iso_incremental[0, :]) *
                                        iso_incremental)


def test_isomap_incremental_adapt():
    # Same as above, only n_max_samples is set to smaller values remaining
    # Manifold is trained adaptively -> must be equal to normal isomap on last samples
    N_per_side = 10
    Npts = N_per_side ** 2
    n_neighbors = 10

    # grid of equidistant points in 2D, n_components = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)), np.float)
    X += np.random.randn(X.shape[0], X.shape[1])

    # select training data
    n_train = 20
    i_train = np.arange(n_train)  # select first n_train samples for training
    i_train.sort()
    i_test = np.setdiff1d(range(Npts), i_train)
    X_train = X[i_train]
    X_test = X[i_test]

    # Incremental Isomap
    clf = manifold.IncrementalIsomap(n_neighbors=n_neighbors, n_components=2,
                                     n_max_samples=2 * n_train,
                                     overflow_mode="adapt")
    clf.fit(X_train)
    for i in range(len(i_test)):
        clf.partial_fit(X_test[i, :].reshape(1, -1))

    iso_incremental = clf.embedding_
    iso_incremental *= 1.0 / np.sqrt(clf.n) * np.sqrt(clf.kernel_pca_.lambdas_)

    # Full isomap only on last samples
    x_readapt = np.concatenate((X[-clf.i - 1:, :], X[-clf.n:-clf.i - 1, :]))
    clf_iso = manifold.Isomap(n_neighbors=n_neighbors, n_components=2)
    iso_full = clf_iso.fit_transform(x_readapt)

    assert_array_almost_equal(iso_full, np.sign(iso_full[0, :]) *  # assure same signs for both embeddings
                                        np.sign(iso_incremental[0, :]) *
                                        iso_incremental)


def test_isomap_reconstruction_error():
    # Same setup as in test_isomap_simple_grid, with an added dimension
    N_per_side = 5
    Npts = N_per_side ** 2
    n_neighbors = Npts - 1

    # grid of equidistant points in 2D, n_components = n_dim
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
            clf = manifold.IncrementalIsomap(n_neighbors=n_neighbors, n_components=2,
                                  eigen_solver=eigen_solver,
                                  path_method=path_method)
            clf.fit(X)

            embedding_unnormed = (clf.embedding_[:clf.n, :] / np.sqrt(clf.n) *
                                    np.sqrt(clf.kernel_pca_.lambdas_))  # convert scale back

            # compute output kernel
            G_iso = neighbors.kneighbors_graph(embedding_unnormed,
                                               n_neighbors,
                                               mode='distance').toarray()

            K_iso = centerer.fit_transform(-0.5 * G_iso ** 2)

            # make sure error agrees
            reconstruction_error = np.linalg.norm(K - K_iso) / Npts
            assert_almost_equal(reconstruction_error,
                                clf.reconstruction_error())


def test_transform():
    n_samples = 200
    n_components = 10
    noise_scale = 0.01

    # Create S-curve dataset
    X, y = datasets.samples_generator.make_s_curve(n_samples)

    # Compute isomap embedding
    iso = manifold.IncrementalIsomap(n_components, 2)
    X_iso = iso.fit_transform(X)

    # Re-embed a noisy version of the points
    rng = np.random.RandomState(0)
    noise = noise_scale * rng.randn(*X.shape)
    X_iso2 = iso.transform(X + noise)

    # Make sure the rms error on re-embedding is comparable to noise_scale
    assert_less(np.sqrt(np.mean((X_iso - X_iso2) ** 2)), 2 * noise_scale)


def test_pipeline():
    # check that Isomap works fine as a transformer in a Pipeline
    # only checks that no error is raised.
    # TODO check that it actually does something useful
    X, y = datasets.make_blobs(random_state=0)
    clf = pipeline.Pipeline(
        [('isomap', manifold.IncrementalIsomap()),
         ('clf', neighbors.KNeighborsClassifier())])
    clf.fit(X, y)
    assert_less(.9, clf.score(X, y))


if __name__ == '__main__':
    import nose
    nose.runmodule()
