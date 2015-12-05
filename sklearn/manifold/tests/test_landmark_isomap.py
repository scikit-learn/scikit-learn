from itertools import product

import numpy as np

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from sklearn import datasets
from sklearn import manifold
from sklearn import neighbors
from sklearn import pipeline
from sklearn import preprocessing
from sklearn.utils import testing

eigen_solvers = ['auto', 'dense', 'arpack']
path_methods = ['auto', 'FW', 'D']


def test_landmarks():
    n_samples = 100
    n_landmarks = n_samples // 3

    X, y = datasets.make_s_curve(n_samples, random_state=0)

    # Should always have a minimum of n_components + 1 landmarks.
    i = manifold.LandmarkIsomap().fit(X)
    testing.assert_less(i.n_components, len(i.landmarks_))

    # Should find n random landmarks.
    landmarks = manifold.LandmarkIsomap(n_landmarks=n_landmarks) \
        .fit(X).landmarks_
    testing.assert_equal(len(landmarks), n_landmarks)

    # Should keep the landmarks passed.
    n_landmarks = 23
    expected_landmarks = np.random.randint(X.shape[0], size=n_landmarks)

    i = manifold.LandmarkIsomap(landmarks=expected_landmarks).fit(X)
    actual_landmarks = i.landmarks_

    expected_landmarks = np.sort(expected_landmarks)
    actual_landmarks = np.sort(actual_landmarks)
    testing.assert_array_equal(actual_landmarks, expected_landmarks)


def test_transform():
    n_samples = 200
    n_landmarks = 20
    n_components = 10
    noise_scale = 0.01

    # Create S-curve dataset
    X, y = datasets.samples_generator.make_s_curve(n_samples, random_state=0)

    # Compute l-isomap embedding
    iso = manifold.LandmarkIsomap(n_components, 2, n_landmarks=n_landmarks)
    X_iso = iso.fit_transform(X)

    # Re-embed a noisy version of the points
    rng = np.random.RandomState(0)
    noise = noise_scale * rng.randn(*X.shape)
    X_iso2 = iso.transform(X + noise)

    # Make sure the rms error on re-embedding is comparable to noise_scale
    testing.assert_less(
        np.sqrt(np.mean((X_iso - X_iso2) ** 2)), 2 * noise_scale)


def test_pipeline():
    # check that Isomap works fine as a transformer in a Pipeline
    # only checks that no error is raised.
    # TODO check that it actually does something useful
    landmarks = list(range(0, 100, 4))
    X, y = datasets.make_blobs(random_state=0)
    clf = pipeline.Pipeline(
        [('l-isomap', manifold.LandmarkIsomap(landmarks=landmarks)),
         ('clf', neighbors.KNeighborsClassifier())])
    clf.fit(X, y)
    testing.assert_less(.9, clf.score(X, y))
