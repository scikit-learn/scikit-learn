import numpy as np
from ..sparsepca import SparsePCA
from numpy.testing import assert_array_almost_equal, assert_equal

def generate_toy_data(n_atoms, n_samples, image_size):
    n_features = image_size[0] * image_size[1]

    np.random.seed(0)
    U = np.random.randn(n_samples, n_atoms)
    V = np.random.randn(n_atoms, n_features)

    centers = [(3, 3), (6, 7), (8, 1)]
    sz = [1, 2, 1]
    for k in range(n_atoms):
        img = np.zeros(image_size)
        xmin, xmax = centers[k][0] - sz[k], centers[k][0] + sz[k]
        ymin, ymax = centers[k][1] - sz[k], centers[k][1] + sz[k]
        img[xmin:xmax][:, ymin:ymax] = 1.0
        V[k, :] = img.ravel()

    # Y is defined by : Y = UV + noise
    Y = np.dot(U, V)
    Y += 0.1 * np.random.randn(Y.shape[0], Y.shape[1])  # Add noise
    return Y, U, V
   
def test_correct_shapes():
    np.random.seed(0)
    X = np.random.randn(12, 10)
    pca = SparsePCA(n_components=8)
    U = pca.fit_transform(X)
    assert_equal(pca.components_.shape, (8, 10))
    assert_equal(U.shape, (12, 8))

def test_fit_transform():
    Y, _, _ = generate_toy_data(3, 10, (8, 8))
    U1 = SparsePCA(n_components=3).fit_transform(Y)
    U2 = SparsePCA(n_components=3).fit(Y).transform(Y)
    assert_array_almost_equal(U1, U2)