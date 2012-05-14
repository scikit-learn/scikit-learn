import numpy as np
import cProfile, pstats

from ..classes import MLPClassifier
from ... import datasets, preprocessing
from numpy.testing import assert_array_almost_equal

def test_xor():
    """Test if MLP is able to learn XOR"""
    X = np.asarray([[0, 0], [0, 1], [1, 0], [1, 1]])
    Y = np.asarray([[0], [1], [1], [0]])

    mlp = MLPClassifier(n_hidden=10, lr=0.1, batch_size=4)

    mlp.fit(X, Y, max_epochs=1000)
    p = mlp.predict(X)

    assert_array_almost_equal(p, Y)

def test_numbers():
    digits = datasets.load_digits()
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))

    data = preprocessing.scale(data)

    classifier = MLPClassifier(n_hidden=10, lr=0.3, batch_size=100)
    classifier.fit(data[:n_samples / 2], digits.target[:n_samples / 2],
        max_epochs=1000)
    expected = digits.target[n_samples / 2:]

    predicted = classifier.predict(data[n_samples / 2:])

    assert_array_almost_equal(predicted, expected)