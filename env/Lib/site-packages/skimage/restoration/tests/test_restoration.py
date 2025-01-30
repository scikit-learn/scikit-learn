import numpy as np
import pytest
from scipy import ndimage as ndi
from scipy.signal import convolve2d, convolve

from skimage import restoration, util
from skimage._shared import filters
from skimage._shared.testing import fetch
from skimage._shared.utils import _supported_float_type
from skimage.color import rgb2gray
from skimage.data import astronaut, camera
from skimage.restoration import uft


test_img = util.img_as_float(camera())


def _get_rtol_atol(dtype):
    rtol = 1e-3
    atol = 0
    if dtype == np.float16:
        rtol = 1e-2
        atol = 1e-3
    elif dtype == np.float32:
        atol = 1e-5
    return rtol, atol


@pytest.mark.parametrize('dtype', [np.float16, np.float32, np.float64])
@pytest.mark.parametrize('ndim', [1, 2, 3])
def test_wiener(dtype, ndim):
    """
    currently only performs pixelwise comparison to
    precomputed result in 2d case.
    """

    rng = np.random.RandomState(0)
    psf = np.ones([5] * ndim, dtype=dtype) / 5**ndim

    # for ndim == 2 use camera (to compare to presaved result)
    if ndim != 2:
        test_img = rng.randint(0, 100, [50] * ndim)
    else:
        test_img = util.img_as_float(camera())

    data = convolve(test_img, psf, 'same')
    data += 0.1 * data.std() * rng.standard_normal(data.shape)
    data = data.astype(dtype, copy=False)
    deconvolved = restoration.wiener(data, psf, 0.05)
    assert deconvolved.dtype == _supported_float_type(dtype)

    if ndim == 2:
        rtol, atol = _get_rtol_atol(dtype)
        path = fetch('restoration/tests/camera_wiener.npy')
        np.testing.assert_allclose(deconvolved, np.load(path), rtol=rtol, atol=atol)

    _, laplacian = uft.laplacian(ndim, data.shape)
    otf = uft.ir2tf(psf, data.shape, is_real=False)
    assert otf.real.dtype == _supported_float_type(dtype)
    deconvolved = restoration.wiener(data, otf, 0.05, reg=laplacian, is_real=False)
    assert deconvolved.real.dtype == _supported_float_type(dtype)
    if ndim == 2:
        np.testing.assert_allclose(
            np.real(deconvolved), np.load(path), rtol=rtol, atol=atol
        )


@pytest.mark.parametrize('dtype', [np.float16, np.float32, np.float64])
def test_unsupervised_wiener(dtype):
    psf = np.ones((5, 5), dtype=dtype) / 25
    data = convolve2d(test_img, psf, 'same')
    seed = 16829302
    # keep old-style RandomState here for compatibility with previously stored
    # reference data in camera_unsup.npy and camera_unsup2.npy
    rng = np.random.RandomState(seed)
    data += 0.1 * data.std() * rng.standard_normal(data.shape)
    data = data.astype(dtype, copy=False)
    deconvolved, _ = restoration.unsupervised_wiener(data, psf, rng=seed)
    restoration.unsupervised_wiener(data, psf, rng=seed)
    float_type = _supported_float_type(dtype)
    assert deconvolved.dtype == float_type

    rtol, atol = _get_rtol_atol(dtype)
    path = fetch('restoration/tests/camera_unsup.npy')
    np.testing.assert_allclose(deconvolved, np.load(path), rtol=rtol, atol=atol)

    _, laplacian = uft.laplacian(2, data.shape)
    otf = uft.ir2tf(psf, data.shape, is_real=False)
    assert otf.real.dtype == _supported_float_type(dtype)
    deconvolved2 = restoration.unsupervised_wiener(
        data,
        otf,
        reg=laplacian,
        is_real=False,
        user_params={
            "callback": lambda x: None,
            "max_num_iter": 200,
            "min_num_iter": 30,
        },
        rng=seed,
    )[0]
    assert deconvolved2.real.dtype == float_type
    path = fetch('restoration/tests/camera_unsup2.npy')
    np.testing.assert_allclose(
        np.real(deconvolved2), np.load(path), rtol=rtol, atol=atol
    )


def test_unsupervised_wiener_deprecated_user_param():
    psf = np.ones((5, 5), dtype=float) / 25
    data = convolve2d(test_img, psf, 'same')
    otf = uft.ir2tf(psf, data.shape, is_real=False)
    _, laplacian = uft.laplacian(2, data.shape)
    restoration.unsupervised_wiener(
        data,
        otf,
        reg=laplacian,
        is_real=False,
        user_params={"max_num_iter": 300, "min_num_iter": 30},
        rng=5,
    )


def test_image_shape():
    """Test that shape of output image in deconvolution is same as input.

    This addresses issue #1172.
    """
    point = np.zeros((5, 5), float)
    point[2, 2] = 1.0
    psf = filters.gaussian(point, sigma=1.0, mode='reflect')
    # image shape: (45, 45), as reported in #1172
    image = util.img_as_float(camera()[65:165, 215:315])  # just the face
    image_conv = ndi.convolve(image, psf)
    deconv_sup = restoration.wiener(image_conv, psf, 1)
    deconv_un = restoration.unsupervised_wiener(image_conv, psf)[0]
    # test the shape
    np.testing.assert_equal(image.shape, deconv_sup.shape)
    np.testing.assert_equal(image.shape, deconv_un.shape)
    # test the reconstruction error
    sup_relative_error = np.abs(deconv_sup - image) / image
    un_relative_error = np.abs(deconv_un - image) / image
    np.testing.assert_array_less(np.median(sup_relative_error), 0.1)
    np.testing.assert_array_less(np.median(un_relative_error), 0.1)


@pytest.mark.parametrize('ndim', [1, 2, 3])
def test_richardson_lucy(ndim):
    psf = np.ones([5] * ndim, dtype=float) / 5**ndim
    if ndim != 2:
        test_img = np.random.randint(0, 100, [30] * ndim)
    else:
        test_img = util.img_as_float(camera())
    data = convolve(test_img, psf, 'same')

    rng = np.random.RandomState(0)
    data += 0.1 * data.std() * rng.standard_normal(data.shape)
    deconvolved = restoration.richardson_lucy(data, psf, num_iter=5)

    if ndim == 2:
        path = fetch('restoration/tests/camera_rl.npy')
        np.testing.assert_allclose(deconvolved, np.load(path), rtol=1e-3)


@pytest.mark.parametrize('dtype_image', [np.float16, np.float32, np.float64])
@pytest.mark.parametrize('dtype_psf', [np.float32, np.float64])
def test_richardson_lucy_filtered(dtype_image, dtype_psf):
    if dtype_image == np.float64:
        atol = 1e-8
    else:
        atol = 1e-5
    test_img_astro = rgb2gray(astronaut())

    psf = np.ones((5, 5), dtype=dtype_psf) / 25
    data = convolve2d(test_img_astro, psf, 'same')
    data = data.astype(dtype_image, copy=False)

    deconvolved = restoration.richardson_lucy(data, psf, 5, filter_epsilon=1e-6)
    assert deconvolved.dtype == _supported_float_type(data.dtype)

    path = fetch('restoration/tests/astronaut_rl.npy')
    np.testing.assert_allclose(deconvolved, np.load(path), rtol=1e-3, atol=atol)
