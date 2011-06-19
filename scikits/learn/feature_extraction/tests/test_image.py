# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np
import scipy as sp
from scipy import ndimage

from nose.tools import assert_equal

from ..image import img_to_graph, grid_to_graph
from ..image import extract_patches_2d, reconstruct_patches
from ...utils.graph import cs_graph_components

def test_img_to_graph():
    x, y = np.mgrid[:4, :4] - 10
    grad_x = img_to_graph(x)
    grad_y = img_to_graph(y)
    assert_equal(grad_x.nnz, grad_y.nnz)
    # Negative elements are the diagonal: the elements of the original
    # image. Positive elements are the values of the gradient, they
    # should all be equal on grad_x and grad_y
    np.testing.assert_array_equal(grad_x.data[grad_x.data > 0],
                                  grad_y.data[grad_y.data > 0])


def test_connect_regions():
    lena = sp.lena()
    for thr in (50, 150):
        mask = lena > thr
        graph = img_to_graph(lena, mask)
        assert_equal(ndimage.label(mask)[1], cs_graph_components(graph)[0])


def test_connect_regions_with_grid():
    lena = sp.lena()
    mask = lena > 50
    graph = grid_to_graph(*lena.shape, **{'mask' : mask})
    assert_equal(ndimage.label(mask)[1], cs_graph_components(graph)[0])

    mask = lena > 150
    graph = grid_to_graph(*lena.shape, **{'mask' : mask, 'dtype' : None})
    assert_equal(ndimage.label(mask)[1], cs_graph_components(graph)[0])


def _downsampled_lena():
    lena = sp.lena()
    lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]
    lena /= 4.0
    return lena

def _orange_lena():
    lena = _downsampled_lena()
    lena_color = np.zeros(lena.shape + (3,))
    lena_color[:, :, 0] = 256 - lena
    lena_color[:, :, 1] = 256 - lena / 2
    lena_color[:, :, 2] = 256 - lena / 4
    return lena_color

def _make_images():
    lena = _downsampled_lena()

    # make a collection of lenas
    images = np.zeros((3,) + lena.shape)
    images[0] = lena
    images[1] = lena + 1
    images[2] = lena + 2
    return images


def test_extract_patches_all():
    lena = _downsampled_lena()
    i_h, i_w = lena.shape
    p_h, p_w = 16, 16
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)

    patches = extract_patches_2d(lena, (i_h, i_w), (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))


def test_extract_patches_all_color():
    lena = _orange_lena()
    i_h, i_w = lena.shape[:2]
    p_h, p_w = 16, 16
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)

    patches = extract_patches_2d(lena, (i_h, i_w, 3), (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w, 3))


def test_extract_patches_all_rect():
    lena = _downsampled_lena()
    lena = lena[:, 32:97]
    i_h, i_w = lena.shape
    p_h, p_w = 16, 12    
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)

    patches = extract_patches_2d(lena, (i_h, i_w), (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))


def test_extract_patches_max_patches():
    lena = _downsampled_lena()
    i_h, i_w = lena.shape
    p_h, p_w = 16, 16

    patches = extract_patches_2d(lena, (i_h, i_w), (p_h, p_w), max_patches=100)
    assert_equal(patches.shape, (100, p_h, p_w))


def test_reconstruct_patches_perfect():
    lena = _downsampled_lena()
    i_h, i_w = lena.shape
    p_h, p_w = 16, 16

    patches = extract_patches_2d(lena, (i_h, i_w), (p_h, p_w))
    lena_reconstructed = reconstruct_patches(patches, (i_h, i_w), (p_h, p_w))
    np.testing.assert_array_equal(lena, lena_reconstructed)

def test_reconstruct_patches_perfect_color():
    lena = _orange_lena()
    p_h, p_w = 16, 16

    patches = extract_patches_2d(lena, lena.shape, (p_h, p_w))
    lena_reconstructed = reconstruct_patches(patches, lena.shape, (p_h, p_w))
    np.testing.assert_array_equal(lena, lena_reconstructed)
