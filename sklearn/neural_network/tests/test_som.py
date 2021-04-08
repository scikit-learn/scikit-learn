"""Testing for K-means"""

import numpy as np
from scipy import sparse as sp

import pytest

from sklearn.utils._testing import assert_almost_equal
from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.neural_network import SOM


@pytest.mark.parametrize("topology", ["rectangular", "hexagonal"])
def test_som_setup(topology):
    som = SOM(8, 8, topology=topology)
    som._randomize_weights(2) # 2 features
    for i in range(8):
        for j in range(8):
            assert (som._weights[i, j] < 1.0).all()

@pytest.mark.parametrize("topology", ["rectangular", "hexagonal"])
def test_euclidean_distance(topology):
    som = SOM(8, 8, topology=topology)
    x = np.zeros((1, 2))
    w = np.ones((2, 2, 2))
    d = som._euclidean(x, w)
    assert_array_almost_equal(d, [[1.41421356, 1.41421356],
                                  [1.41421356, 1.41421356]])

@pytest.mark.parametrize("topology", ["rectangular", "hexagonal"])
def test_cosine_distance(topology):
    som = SOM(8, 8, topology=topology)
    x = np.zeros((1, 2))
    w = np.ones((2, 2, 2))
    d = som._cosine(x, w)
    assert_array_almost_equal(d, [[1., 1.],
                                  [1., 1.]])

@pytest.mark.parametrize("topology", ["rectangular", "hexagonal"])
def test_manhattan_distance(topology):
    som = SOM(8, 8, topology=topology)
    x = np.zeros((1, 2))
    w = np.ones((2, 2, 2))
    d = som._manhattan(x, w)
    assert_array_almost_equal(d, [[2., 2.],
                                  [2., 2.]])

@pytest.mark.parametrize("topology", ["rectangular", "hexagonal"])
def test_chebyshev_distance(topology):
    som = SOM(8, 8, topology=topology)
    x = np.array([1, 3])
    w = np.ones((2, 2, 2))
    d = som._chebyshev(x, w)
    assert_array_almost_equal(d, [[2., 2.],
                                  [2., 2.]])

@pytest.mark.parametrize("topology", ["rectangular", "hexagonal"])
def test_som_predict(topology):
    som = SOM(2, 2, topology=topology)

    y1 = np.array([[0.1, 0.4]])
    y2 = np.array([[0.1, 7.8]])
    y3 = np.array([[7.6, 0.4]])
    y4 = np.array([[7.6, 7.8]])

    som.set_weights(np.array([[[0.25, 0.25], [0.25, 7.75]],
                             [[7.75, 0.25], [7.75, 7.75]]]))

    # predict method requires fitting
    q1 = som.quantization(y1)[-1]
    q2 = som.quantization(y2)[-1]
    q3 = som.quantization(y3)[-1]
    q4 = som.quantization(y4)[-1]

    # first group
    assert_array_equal(q1, [0.25, 0.25])
    # second group
    assert_array_equal(q2, [0.25, 7.75])
    # third group
    assert_array_equal(q3, [7.75, 0.25])
    # fourth group
    assert_array_equal(q4, [7.75, 7.75])

@pytest.mark.parametrize("topology", ["rectangular", "hexagonal"])
@pytest.mark.parametrize("neighborhood_function", ["gaussian", "mexican_hat", "bubble", "triangle"])
@pytest.mark.parametrize("activation_distance", ["euclidean", "cosine", "manhattan", "chebyshev"])
@pytest.mark.parametrize("state", [1, 10, 1048, 1000000])
def test_som_random_seed(topology, neighborhood_function, activation_distance, state):
    som1 = SOM(4, 4, sigma=1.0, learning_rate=0.5, max_iter=100,
               topology=topology, neighborhood_function=neighborhood_function,
               activation_distance=activation_distance, random_state=state)
    som2 = SOM(4, 4, sigma=1.0, learning_rate=0.5, max_iter=100,
               topology=topology, neighborhood_function=neighborhood_function,
               activation_distance=activation_distance, random_state=state)

    som1._randomize_weights(2) # 2 features
    som2._randomize_weights(2) # 2 features

    # same initialization
    assert_array_almost_equal(som1._weights, som2._weights)

    X = np.random.rand(100, 2)

    som1.fit(X)
    som2.fit(X)

    # same result
    assert_array_almost_equal(som1._weights, som2._weights)
