import numpy as np
import scipy as sp
import pytest

from skimage.graph._graph import pixel_graph, central_pixel

mask = np.array([[1, 0, 0], [0, 1, 1], [0, 1, 0]], dtype=bool)
image = np.random.default_rng().random(mask.shape)


def test_small_graph():
    g, n = pixel_graph(mask, connectivity=2)
    assert g.shape == (4, 4)
    assert len(g.data) == 8
    np.testing.assert_allclose(np.unique(g.data), [1, np.sqrt(2)])
    np.testing.assert_array_equal(n, [0, 4, 5, 7])


def test_pixel_graph_return_type():
    g, n = pixel_graph(mask, connectivity=2)
    assert isinstance(g, sp.sparse.csr_matrix)

    g, n = pixel_graph(mask, connectivity=2, sparse_type="matrix")
    assert isinstance(g, sp.sparse.csr_matrix)

    g, n = pixel_graph(mask, connectivity=2, sparse_type="array")
    assert isinstance(g, sp.sparse.csr_array)

    with pytest.raises(ValueError, match="`sparse_type` must be 'array' or 'matrix'"):
        pixel_graph(mask, connectivity=2, sparse_type="unknown")


@pytest.mark.parametrize("sparse_type", ["matrix", "array"])
def test_central_pixel(sparse_type):
    g, n = pixel_graph(mask, connectivity=2, sparse_type=sparse_type)
    px, ds = central_pixel(g, n, shape=mask.shape)
    np.testing.assert_array_equal(px, (1, 1))
    s2 = np.sqrt(2)
    np.testing.assert_allclose(ds, [s2 * 3 + 2, s2 + 2, s2 * 2 + 2, s2 * 2 + 2])

    # test raveled coordinate
    px, _ = central_pixel(g, n)
    assert px == 4

    # test no nodes given
    px, _ = central_pixel(g)
    assert px == 1


@pytest.mark.parametrize("sparse_type", ["matrix", "array"])
def test_edge_function(sparse_type):
    def edge_func(values_src, values_dst, distances):
        return np.abs(values_src - values_dst) + distances

    g, n = pixel_graph(
        image,
        mask=mask,
        connectivity=2,
        edge_function=edge_func,
        sparse_type=sparse_type,
    )
    s2 = np.sqrt(2)
    np.testing.assert_allclose(g[0, 1], np.abs(image[0, 0] - image[1, 1]) + s2)
    np.testing.assert_allclose(g[1, 2], np.abs(image[1, 1] - image[1, 2]) + 1)
    np.testing.assert_array_equal(n, [0, 4, 5, 7])


@pytest.mark.parametrize("sparse_type", ["matrix", "array"])
def test_default_edge_func(sparse_type):
    g, n = pixel_graph(image, spacing=np.array([0.78, 0.78]), sparse_type=sparse_type)
    num_edges = len(g.data) // 2  # each edge appears in both directions
    assert num_edges == 12  # lattice in a (3, 3) grid
    np.testing.assert_almost_equal(g[0, 1], 0.78 * np.abs(image[0, 0] - image[0, 1]))
    np.testing.assert_array_equal(n, np.arange(image.size))


@pytest.mark.parametrize("sparse_type", ["matrix", "array"])
def test_no_mask_with_edge_func(sparse_type):
    """Ensure function `pixel_graph` runs when passing `edge_function` but not `mask`."""
    image = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    def func(x, y, z):
        return np.abs(x - y) * 0.5

    expected_g = (
        np.array(
            [
                [0.0, 1.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 1.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 3.0, 0.0, 0.0],
                [0.0, 3.0, 0.0, 1.0, 0.0, 1.0, 0.0, 3.0, 0.0],
                [0.0, 0.0, 3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 3.0],
                [0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 1.0, 0.0, 1.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 1.0, 0.0],
            ]
        )
        * 0.5
    )

    g, n = pixel_graph(image, edge_function=func, sparse_type=sparse_type)
    np.testing.assert_array_equal(n, np.arange(image.size))
    np.testing.assert_array_equal(g.toarray(), expected_g)
