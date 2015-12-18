import numpy as np

from sklearn import datasets
from sklearn import manifold
from sklearn import neighbors
from sklearn import pipeline
from sklearn.utils.testing import (assert_less, assert_equal, assert_true,
                                   assert_raises)


def test_landmarks():
    X, y = datasets.make_blobs(random_state=0)

    # Should have at least one landmark.
    for method in ('random', 'min-max'):
        isomap = manifold.Isomap(landmarks='auto',
                                 landmark_selection_method=method).fit(X, y)
        assert_less(0, len(isomap.landmarks_))

    expected_n_landmarks = 7
    isomap = manifold.Isomap(landmarks=expected_n_landmarks,
                             landmark_selection_method='random').fit(X, y)
    assert_equal(len(isomap.landmarks_), expected_n_landmarks)

    expected_n_landmarks = 11
    isomap = manifold.Isomap(landmarks=expected_n_landmarks,
                             landmark_selection_method='min-max').fit(X, y)
    assert_equal(len(isomap.landmarks_), expected_n_landmarks)

    expected_landmarks = np.array([0, 2, 4, 10])
    isomap = manifold.Isomap(landmarks=expected_landmarks).fit(X, y)
    assert_true(isomap.landmarks_ is expected_landmarks,
                'L-Isomap did not preserve the landmarks given.')

    # Assert that Isomap raises errors for invalid states.
    iso = manifold.Isomap(landmarks='invalid')
    assert_raises(ValueError, iso.fit, X, y)

    iso = manifold.Isomap(landmarks='auto', landmark_selection_method='invalid')
    assert_raises(ValueError, iso.fit, X, y)


def test_transform():
    n_samples = 200
    n_components = 10
    noise_scale = 0.01

    # Create S-curve dataset
    X, y = datasets.make_s_curve(n_samples, random_state=0)

    # Compute isomap embedding
    iso = manifold.Isomap(n_components, 2, landmarks=n_samples // 10)
    X_iso = iso.fit_transform(X)

    # Re-embed a noisy version of the points
    rng = np.random.RandomState(0)
    noise = noise_scale * rng.randn(*X.shape)
    X_iso2 = iso.transform(X + noise)

    # Make sure the rms error on re-embedding is comparable to noise_scale
    assert_less(np.sqrt(np.mean((X_iso - X_iso2) ** 2)), 2 * noise_scale)


def test_pipeline():
    # check that L-Isomap works fine as a transformer
    # in a Pipeline only checks that no error is raised.
    # TODO check that it actually does something useful
    X, y = datasets.make_blobs(random_state=0)
    clf = pipeline.Pipeline(
            [('isomap', manifold.Isomap(landmarks='auto')),
             ('clf', neighbors.KNeighborsClassifier())])
    clf.fit(X, y)
    assert_less(.9, clf.score(X, y))
