# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD 3 clause

import numpy as np
import scipy as sp
from scipy import ndimage
from scipy.sparse.csgraph import connected_components

from numpy.testing import assert_raises

from sklearn.feature_extraction.image import (
    img_to_graph, grid_to_graph, extract_patches_2d,
    reconstruct_from_patches_2d, PatchExtractor, extract_patches)
from sklearn.utils.testing import SkipTest, assert_equal, assert_true
from sklearn.utils.fixes import sp_version

if sp_version < (0, 12):
    raise SkipTest("Skipping because SciPy version earlier than 0.12.0 and "
                   "thus does not include the scipy.misc.face() image.")


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
    # Checking that the function works with graphs containing no edges
    size = 2
    roi_size = 1
    # Generating two convex parts with one vertex
    # Thus, edges will be empty in _to_graph
    mask = np.zeros((size, size), dtype=np.bool)
    mask[0:roi_size, 0:roi_size] = True
    mask[-roi_size:, -roi_size:] = True
    mask = mask.reshape(size ** 2)
    A = grid_to_graph(n_x=size, n_y=size, mask=mask, return_as=np.ndarray)
    assert_true(connected_components(A)[0] == 2)

    # Checking that the function works whatever the type of mask is
    mask = np.ones((size, size), dtype=np.int16)
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask)
    assert_true(connected_components(A)[0] == 1)

    # Checking dtype of the graph
    mask = np.ones((size, size))
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask, dtype=np.bool)
    assert_true(A.dtype == np.bool)
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask, dtype=np.int)
    assert_true(A.dtype == np.int)
    A = grid_to_graph(n_x=size, n_y=size, n_z=size, mask=mask,
                      dtype=np.float64)
    assert_true(A.dtype == np.float64)


def test_connect_regions():
    try:
        face = sp.face(gray=True)
    except AttributeError:
        # Newer versions of scipy have face in misc
        from scipy import misc
        face = misc.face(gray=True)
    for thr in (50, 150):
        mask = face > thr
        graph = img_to_graph(face, mask)
        assert_equal(ndimage.label(mask)[1], connected_components(graph)[0])


def test_connect_regions_with_grid():
    try:
        face = sp.face(gray=True)
    except AttributeError:
        # Newer versions of scipy have face in misc
        from scipy import misc
        face = misc.face(gray=True)
    mask = face > 50
    graph = grid_to_graph(*face.shape, mask=mask)
    assert_equal(ndimage.label(mask)[1], connected_components(graph)[0])

    mask = face > 150
    graph = grid_to_graph(*face.shape, mask=mask, dtype=None)
    assert_equal(ndimage.label(mask)[1], connected_components(graph)[0])


def _downsampled_face():
    try:
        face = sp.face(gray=True)
    except AttributeError:
        # Newer versions of scipy have face in misc
        from scipy import misc
        face = misc.face(gray=True)
    face = face.astype(np.float32)
    face = (face[::2, ::2] + face[1::2, ::2] + face[::2, 1::2]
            + face[1::2, 1::2])
    face = (face[::2, ::2] + face[1::2, ::2] + face[::2, 1::2]
            + face[1::2, 1::2])
    face = face.astype(np.float32)
    face /= 16.0
    return face


def _orange_face(face=None):
    face = _downsampled_face() if face is None else face
    face_color = np.zeros(face.shape + (3,))
    face_color[:, :, 0] = 256 - face
    face_color[:, :, 1] = 256 - face / 2
    face_color[:, :, 2] = 256 - face / 4
    return face_color


def _make_images(face=None):
    face = _downsampled_face() if face is None else face
    # make a collection of faces
    images = np.zeros((3,) + face.shape)
    images[0] = face
    images[1] = face + 1
    images[2] = face + 2
    return images

downsampled_face = _downsampled_face()
orange_face = _orange_face(downsampled_face)
face_collection = _make_images(downsampled_face)


def test_extract_patches_all():
    face = downsampled_face
    i_h, i_w = face.shape
    p_h, p_w = 16, 16
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)
    patches = extract_patches_2d(face, (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))


def test_extract_patches_all_color():
    face = orange_face
    i_h, i_w = face.shape[:2]
    p_h, p_w = 16, 16
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)
    patches = extract_patches_2d(face, (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w, 3))


def test_extract_patches_all_rect():
    face = downsampled_face
    face = face[:, 32:97]
    i_h, i_w = face.shape
    p_h, p_w = 16, 12
    expected_n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)

    patches = extract_patches_2d(face, (p_h, p_w))
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))


def test_extract_patches_max_patches():
    face = downsampled_face
    i_h, i_w = face.shape
    p_h, p_w = 16, 16

    patches = extract_patches_2d(face, (p_h, p_w), max_patches=100)
    assert_equal(patches.shape, (100, p_h, p_w))

    expected_n_patches = int(0.5 * (i_h - p_h + 1) * (i_w - p_w + 1))
    patches = extract_patches_2d(face, (p_h, p_w), max_patches=0.5)
    assert_equal(patches.shape, (expected_n_patches, p_h, p_w))

    assert_raises(ValueError, extract_patches_2d, face, (p_h, p_w),
                  max_patches=2.0)
    assert_raises(ValueError, extract_patches_2d, face, (p_h, p_w),
                  max_patches=-1.0)


def test_reconstruct_patches_perfect():
    face = downsampled_face
    p_h, p_w = 16, 16

    patches = extract_patches_2d(face, (p_h, p_w))
    face_reconstructed = reconstruct_from_patches_2d(patches, face.shape)
    np.testing.assert_array_almost_equal(face, face_reconstructed)


def test_reconstruct_patches_perfect_color():
    face = orange_face
    p_h, p_w = 16, 16

    patches = extract_patches_2d(face, (p_h, p_w))
    face_reconstructed = reconstruct_from_patches_2d(patches, face.shape)
    np.testing.assert_array_almost_equal(face, face_reconstructed)


def test_patch_extractor_fit():
    faces = face_collection
    extr = PatchExtractor(patch_size=(8, 8), max_patches=100, random_state=0)
    assert_true(extr == extr.fit(faces))


def test_patch_extractor_max_patches():
    faces = face_collection
    i_h, i_w = faces.shape[1:3]
    p_h, p_w = 8, 8

    max_patches = 100
    expected_n_patches = len(faces) * max_patches
    extr = PatchExtractor(patch_size=(p_h, p_w), max_patches=max_patches,
                          random_state=0)
    patches = extr.transform(faces)
    assert_true(patches.shape == (expected_n_patches, p_h, p_w))

    max_patches = 0.5
    expected_n_patches = len(faces) * int((i_h - p_h + 1) * (i_w - p_w + 1)
                                          * max_patches)
    extr = PatchExtractor(patch_size=(p_h, p_w), max_patches=max_patches,
                          random_state=0)
    patches = extr.transform(faces)
    assert_true(patches.shape == (expected_n_patches, p_h, p_w))


def test_patch_extractor_max_patches_default():
    faces = face_collection
    extr = PatchExtractor(max_patches=100, random_state=0)
    patches = extr.transform(faces)
    assert_equal(patches.shape, (len(faces) * 100, 19, 25))


def test_patch_extractor_all_patches():
    faces = face_collection
    i_h, i_w = faces.shape[1:3]
    p_h, p_w = 8, 8
    expected_n_patches = len(faces) * (i_h - p_h + 1) * (i_w - p_w + 1)
    extr = PatchExtractor(patch_size=(p_h, p_w), random_state=0)
    patches = extr.transform(faces)
    assert_true(patches.shape == (expected_n_patches, p_h, p_w))


def test_patch_extractor_color():
    faces = _make_images(orange_face)
    i_h, i_w = faces.shape[1:3]
    p_h, p_w = 8, 8
    expected_n_patches = len(faces) * (i_h - p_h + 1) * (i_w - p_w + 1)
    extr = PatchExtractor(patch_size=(p_h, p_w), random_state=0)
    patches = extr.transform(faces)
    assert_true(patches.shape == (expected_n_patches, p_h, p_w, 3))


def test_extract_patches_strided():

    image_shapes_1D = [(10,), (10,), (11,), (10,)]
    patch_sizes_1D = [(1,), (2,), (3,), (8,)]
    patch_steps_1D = [(1,), (1,), (4,), (2,)]

    expected_views_1D = [(10,), (9,), (3,), (2,)]
    last_patch_1D = [(10,), (8,), (8,), (2,)]

    image_shapes_2D = [(10, 20), (10, 20), (10, 20), (11, 20)]
    patch_sizes_2D = [(2, 2), (10, 10), (10, 11), (6, 6)]
    patch_steps_2D = [(5, 5), (3, 10), (3, 4), (4, 2)]

    expected_views_2D = [(2, 4), (1, 2), (1, 3), (2, 8)]
    last_patch_2D = [(5, 15), (0, 10), (0, 8), (4, 14)]

    image_shapes_3D = [(5, 4, 3), (3, 3, 3), (7, 8, 9), (7, 8, 9)]
    patch_sizes_3D = [(2, 2, 3), (2, 2, 2), (1, 7, 3), (1, 3, 3)]
    patch_steps_3D = [(1, 2, 10), (1, 1, 1), (2, 1, 3), (3, 3, 4)]

    expected_views_3D = [(4, 2, 1), (2, 2, 2), (4, 2, 3), (3, 2, 2)]
    last_patch_3D = [(3, 2, 0), (1, 1, 1), (6, 1, 6), (6, 3, 4)]

    image_shapes = image_shapes_1D + image_shapes_2D + image_shapes_3D
    patch_sizes = patch_sizes_1D + patch_sizes_2D + patch_sizes_3D
    patch_steps = patch_steps_1D + patch_steps_2D + patch_steps_3D
    expected_views = expected_views_1D + expected_views_2D + expected_views_3D
    last_patches = last_patch_1D + last_patch_2D + last_patch_3D

    for (image_shape, patch_size, patch_step, expected_view,
         last_patch) in zip(image_shapes, patch_sizes, patch_steps,
                            expected_views, last_patches):
        image = np.arange(np.prod(image_shape)).reshape(image_shape)
        patches = extract_patches(image, patch_shape=patch_size,
                                  extraction_step=patch_step)

        ndim = len(image_shape)

        assert_true(patches.shape[:ndim] == expected_view)
        last_patch_slices = [slice(i, i + j, None) for i, j in
                             zip(last_patch, patch_size)]
        assert_true((patches[[slice(-1, None, None)] * ndim] ==
                    image[last_patch_slices].squeeze()).all())


def test_extract_patches_square():
    # test same patch size for all dimensions
    face = downsampled_face
    i_h, i_w = face.shape
    p = 8
    expected_n_patches = ((i_h - p + 1), (i_w - p + 1))
    patches = extract_patches(face, patch_shape=p)
    assert_true(patches.shape == (expected_n_patches[0], expected_n_patches[1],
                                  p, p))


def test_width_patch():
    # width and height of the patch should be less than the image
    x = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert_raises(ValueError, extract_patches_2d, x, (4, 1))
    assert_raises(ValueError, extract_patches_2d, x, (1, 4))
