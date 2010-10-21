import numpy as np
import numpy.linalg as la

from numpy.testing import assert_array_almost_equal, assert_array_equal, \
                          assert_almost_equal, assert_equal

from scikits.learn.preprocessing import Scaler, scale, Normalizer, \
                                        LengthNormalizer, Binarizer

np.random.seed(0)

def test_scaler():
    """Test scaling of dataset along all axis
    """
    X = np.random.randn(4, 5)

    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 5*[0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), 5*[1.0])
    # Check that X has not been copied
    assert X_scaled is X

    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_array_almost_equal(X_scaled.mean(axis=0), 5*[0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), 5*[1.0])
    # Check that X has not been copied
    assert X_scaled is not X

    X_scaled = scale(X, axis=1, with_std=False)
    assert_array_almost_equal(X_scaled.mean(axis=1), 4*[0.0])
    X_scaled = scale(X, axis=1, with_std=True)
    assert_array_almost_equal(X_scaled.std(axis=1), 4*[1.0])
    # Check that the data hasn't been modified

def test_normalizer():
    X = np.random.randn(4, 5)

    normalizer = Normalizer()
    X_norm = normalizer.transform(X, copy=False)
    assert_array_almost_equal(X_norm.sum(axis=1), np.ones(X.shape[0]))
    assert X_norm is X

    normalizer = Normalizer()
    X_norm = normalizer.transform(X, copy=True)
    assert_array_almost_equal(X_norm.sum(axis=1), np.ones(X.shape[0]))
    assert X_norm is not X

def test_length_normalizer():
    X = np.random.randn(4, 5)

    normalizer = LengthNormalizer()
    X_norm1 = normalizer.transform(X, copy=True)
    assert X_norm1 is not X

    normalizer = LengthNormalizer()
    X_norm2 = normalizer.transform(X, copy=False)
    assert X_norm2 is X

    for X_norm in (X_norm1, X_norm2):
        for i in xrange(len(X_norm)):
            assert_almost_equal(la.norm(X_norm[i]), 1.0)

def test_binarizer():
    X = np.array([[1, 0, 5],
                  [2, 3, 0]])

    binarizer = Binarizer(threshold=2.0)
    X_bin = binarizer.transform(X, copy=True)
    assert_equal(np.sum(X_bin==0), 4)
    assert_equal(np.sum(X_bin==1), 2)

    binarizer = Binarizer()
    X_bin = binarizer.transform(X, copy=True)
    assert_equal(np.sum(X_bin==0), 2)
    assert_equal(np.sum(X_bin==1), 4)
    assert X_bin is not X

    binarizer = Binarizer()
    X_bin = binarizer.transform(X, copy=False)
    assert_equal(np.sum(X_bin==0), 2)
    assert_equal(np.sum(X_bin==1), 4)
    assert X_bin is X


