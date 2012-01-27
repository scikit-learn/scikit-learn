# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np
import scipy as sp
from scipy import ndimage

from nose.tools import assert_equal
from numpy.testing import assert_raises

from ..image import img_to_graph, grid_to_graph
from ..image import extract_patches_2d, reconstruct_from_patches_2d, \
                    PatchExtractor
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


def test_grid_to_graph():
    #Checking that the function works with graphs containing no edges
    size = 2
    roi_size = 1
    # Generating two convex parts with one vertex
    # Thus, edges will be empty in _to_graph
    mask = np.zeros((size, size), dtype=np.bool)
    mask[0:roi_size, 0:roi_size] = True
    mask[-roi_size:, -roi_size:] = True
    mask = mask.reshape(size ** 2)
    A = grid_to_graph(n_x=size, n_y=size, mask=mask, return_as=np.ndarray)
    assert(cs_graph_components(A)[0] == 2)

    # Checking that the function works whatever the type of mask is
    mask = np.ones((size, size), dtype=np.int16)
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask)
    assert(cs_graph_components(A)[0] == 1)

    # Checking dtype of the graph
    mask = np.ones((size, size))
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask, dtype=np.bool)
    assert A.dtype == np.bool
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask, dtype=np.int)
    assert A.dtype == np.int
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask, dtype=np.float)
    assert A.dtype == np.float


def test_connect_regions():
    lena = sp.misc.lena()
    for thr in (50, 150):
        mask = lena > thr
        graph = img_to_graph(lena, mask)
        assert_equal(ndimage.label(mask)[1], cs_graph_components(graph)[0])


def test_connect_regions_with_grid():
    lena = sp.misc.lena()
    mask = lena > 50
    graph = grid_to_graph(*lena.shape, mask=mask)
    assert_equal(ndimage.label(mask)[1], cs_graph_components(graph)[0])

    mask = lena > 150
    graph = grid_to_graph(*lena.shape, mask=mask, dtype=None)
    assert_equal(ndimage.label(mask)[1], cs_graph_components(graph)[0])


def _downsampled_lena():
    lena = sp.misc.lena().astype(np.float32)
    lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + \
           lena[1::2, 1::2]
    lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + \
           lena[1::2, 1::2]
    lena = lena.astype(np.float)
    lena /= 16.0
    return lena


def _orange_lena(lena=None):
    lena = _downsampled_lena() if lena is None else lena
    lena_color = np.zeros(lena.shape + (3,))
    lena_color[:, :, 0] = 256 - lena
    lena_color[:, :, 1] = 256 - lena / 2
    lena_color[:, :, 2] = 256 - lena / 4
    return lena_color


def _make_images(lena=None):
    lena = _downsampled_lena() if lena is None else lena
    # make a collection of lenas
    images = np.zeros((3,) + lena.shape)
    images[0] = lena
    images[1] = lena + 1
    images[2] = lena + 2
    return images

downsampled_lena = _downsampled_lena()
orange_lena = _orange_lena(downsampled_lena)
lena_collection = _make_images(downsampled_lena)


def test_extract_patches_all():
    lena = downsampled_lena
    i_h, i_w = lena.shape
    p_h, p_w = 16, 16
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)
    patches = extract_patches_2d(lena, (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))


def test_extract_patches_all_color():
    lena = orange_lena
    i_h, i_w = lena.shape[:2]
    p_h, p_w = 16, 16
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)
    patches = extract_patches_2d(lena, (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w, 3))


def test_extract_patches_all_rect():
    lena = downsampled_lena
    lena = lena[:, 32:97]
    i_h, i_w = lena.shape
    p_h, p_w = 16, 12
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)

    patches = extract_patches_2d(lena, (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))


def test_extract_patches_max_patches():
    lena = downsampled_lena
    i_h, i_w = lena.shape
    p_h, p_w = 16, 16

    patches = extract_patches_2d(lena, (p_h, p_w), max_patches=100)
    assert_equal(patches.shape, (100, p_h, p_w))

    expected_n_patches = int(0.5 * (i_h - p_h + 1) * (i_w - p_w + 1))
    patches = extract_patches_2d(lena, (p_h, p_w), max_patches=0.5)
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))

    assert_raises(ValueError, extract_patches_2d, lena,
                                                  (p_h, p_w),
                                                  max_patches=2.0)
    assert_raises(ValueError, extract_patches_2d, lena,
                                                  (p_h, p_w),
                                                  max_patches=-1.0)


def test_reconstruct_patches_perfect():
    lena = downsampled_lena
    p_h, p_w = 16, 16

    patches = extract_patches_2d(lena, (p_h, p_w))
    lena_reconstructed = reconstruct_from_patches_2d(patches, lena.shape)
    np.testing.assert_array_equal(lena, lena_reconstructed)


def test_reconstruct_patches_perfect_color():
    lena = orange_lena
    p_h, p_w = 16, 16

    patches = extract_patches_2d(lena, (p_h, p_w))
    lena_reconstructed = reconstruct_from_patches_2d(patches, lena.shape)
    np.testing.assert_array_equal(lena, lena_reconstructed)


def test_patch_extractor_fit():
    lenas = lena_collection
    extr = PatchExtractor(patch_size=(8, 8), max_patches=100, random_state=0)
    assert extr == extr.fit(lenas)


def test_patch_extractor_max_patches():
    lenas = lena_collection
    extr = PatchExtractor(patch_size=(8, 8), max_patches=100, random_state=0)
    patches = extr.transform(lenas)
    assert patches.shape == (len(lenas) * 100, 8, 8)


def test_patch_extractor_all_patches():
    lenas = lena_collection
    i_h, i_w = lenas.shape[1:3]
    p_h, p_w = 8, 8
    expected_n_patches = len(lenas) * (i_h - p_h + 1) * (i_w - p_w + 1)
    extr = PatchExtractor(patch_size=(p_h, p_w), random_state=0)
    patches = extr.transform(lenas)
    assert patches.shape == (expected_n_patches, p_h, p_w)


def test_patch_extractor_color():
    lenas = _make_images(orange_lena)
    i_h, i_w = lenas.shape[1:3]
    p_h, p_w = 8, 8
    expected_n_patches = len(lenas) * (i_h - p_h + 1) * (i_w - p_w + 1)
    extr = PatchExtractor(patch_size=(p_h, p_w), random_state=0)
    patches = extr.transform(lenas)
    assert patches.shape == (expected_n_patches, p_h, p_w, 3)

if __name__ == '__main__':
    import nose
    nose.runmodule()
