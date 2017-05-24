from itertools import product

import numpy as np
from numpy.testing import (assert_almost_equal, assert_array_almost_equal,
                           assert_equal)

from sklearn import datasets
from sklearn import manifold
from sklearn import neighbors
from sklearn import pipeline
from sklearn import preprocessing
from sklearn.utils.testing import (assert_less, assert_true, assert_raises)

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
            clf = manifold.Isomap(n_neighbors=n_neighbors, n_components=2,
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
            clf = manifold.Isomap(n_neighbors=n_neighbors, n_components=2,
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
    n_samples = 200
    n_components = 10
    noise_scale = 0.01

    # Create S-curve dataset
    X, y = datasets.samples_generator.make_s_curve(n_samples, random_state=0)

    # Compute isomap embedding
    iso = manifold.Isomap(n_components, 2)
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
        [('isomap', manifold.Isomap()),
         ('clf', neighbors.KNeighborsClassifier())])
    clf.fit(X, y)
    assert_less(.9, clf.score(X, y))


def test_isomap_clone_bug():
    # regression test for bug reported in #6062
    model = manifold.Isomap()
    for n_neighbors in [10, 15, 20]:
        model.set_params(n_neighbors=n_neighbors)
        model.fit(np.random.rand(50, 2))
        assert_equal(model.nbrs_.n_neighbors,
                     n_neighbors)


def test_landmarks():
    X, y = datasets.make_blobs(random_state=0)

    # Should have at least one landmark.
    for method in ('random', 'min-max'):
        isomap = manifold.Isomap(landmarks='auto',
                                 landmarks_method=method).fit(X, y)
        assert_less(0, len(isomap.landmarks_))

    expected_n_landmarks = 7
    isomap = manifold.Isomap(landmarks=expected_n_landmarks,
                             landmarks_method='random').fit(X, y)
    assert_equal(len(isomap.landmarks_), expected_n_landmarks)

    expected_n_landmarks = 11
    isomap = manifold.Isomap(landmarks=expected_n_landmarks,
                             landmarks_method='min-max').fit(X, y)
    assert_equal(len(isomap.landmarks_), expected_n_landmarks)

    expected_landmarks = np.array([0, 2, 4, 10])
    isomap = manifold.Isomap(landmarks=expected_landmarks).fit(X, y)
    assert_true(isomap.landmarks_ is expected_landmarks,
                'L-Isomap did not preserve the landmarks given.')

    # Assert that Isomap raises errors for invalid states.
    isomap = manifold.Isomap(landmarks='invalid')
    assert_raises(ValueError, isomap.fit, X, y)

    isomap = manifold.Isomap(landmarks='auto', landmarks_method='invalid')
    assert_raises(ValueError, isomap.fit, X, y)


def test_transform_using_landmarks():
    n_samples = 200
    n_components = 10
    noise_scale = 0.01

    # Create S-curve dataset
    X, y = datasets.make_s_curve(n_samples, random_state=0)

    # Compute isomap embedding
    isomap = manifold.Isomap(n_components, 2, landmarks=n_samples // 10)
    X_iso = isomap.fit_transform(X)

    # Re-embed a noisy version of the points
    rng = np.random.RandomState(0)
    noise = noise_scale * rng.randn(*X.shape)
    X_iso2 = isomap.transform(X + noise)

    # Make sure the rms error on re-embedding is comparable to noise_scale
    assert_less(np.sqrt(np.mean((X_iso - X_iso2) ** 2)), 2 * noise_scale)


def test_pipeline_using_landmarks():
    # check that L-Isomap works fine as a transformer
    # in a Pipeline only checks that no error is raised.
    # TODO check that it actually does something useful
    X, y = datasets.make_blobs(random_state=0)
    clf = pipeline.Pipeline(
        [('isomap', manifold.Isomap(landmarks='auto')),
         ('clf', neighbors.KNeighborsClassifier())])
    clf.fit(X, y)
    assert_less(.9, clf.score(X, y))
