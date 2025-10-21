import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.datasets import load_iris
from sklearn.decomposition import PCA
from sklearn.manifold import ClassicalMDS
from sklearn.metrics import euclidean_distances


def test_classical_mds_equivalent_to_pca():
    X, _ = load_iris(return_X_y=True)

    cmds = ClassicalMDS(n_components=2, metric="euclidean")
    pca = PCA(n_components=2)

    Z1 = cmds.fit_transform(X)
    Z2 = pca.fit_transform(X)

    # Swap the signs if necessary
    for comp in range(2):
        if Z1[0, comp] < 0 and Z2[0, comp] > 0:
            Z2[:, comp] *= -1

    assert_allclose(Z1, Z2)

    assert_allclose(np.sqrt(cmds.eigenvalues_), pca.singular_values_)


def test_classical_mds_equivalent_on_data_and_distances():
    X, _ = load_iris(return_X_y=True)

    cmds = ClassicalMDS(n_components=2, metric="euclidean")
    Z1 = cmds.fit_transform(X)

    cmds = ClassicalMDS(n_components=2, metric="precomputed")
    Z2 = cmds.fit_transform(euclidean_distances(X))

    assert_allclose(Z1, Z2)


def test_classical_mds_wrong_inputs():
    # Non-symmetric input
    dissim = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
    with pytest.raises(ValueError, match="Array must be symmetric"):
        ClassicalMDS(metric="precomputed").fit(dissim)

    # Non-square input
    dissim = np.array([[0, 1, 2], [3, 4, 5]])
    with pytest.raises(ValueError, match="array must be 2-dimensional and square"):
        ClassicalMDS(metric="precomputed").fit(dissim)


def test_classical_mds_metric_params():
    X, _ = load_iris(return_X_y=True)

    cmds = ClassicalMDS(n_components=2, metric="euclidean")
    Z1 = cmds.fit_transform(X)

    cmds = ClassicalMDS(n_components=2, metric="minkowski", metric_params={"p": 2})
    Z2 = cmds.fit_transform(X)

    assert_allclose(Z1, Z2)

    cmds = ClassicalMDS(n_components=2, metric="minkowski", metric_params={"p": 1})
    Z3 = cmds.fit_transform(X)

    assert not np.allclose(Z1, Z3)
