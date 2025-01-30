import numpy as np
import pytest

from skimage.io import imread, imsave, use_plugin, reset_plugins, plugin_order
from skimage._shared import testing


pytest.importorskip('SimpleITK')


@pytest.fixture(autouse=True)
def use_simpleitk_plugin():
    """Ensure that SimpleITK plugin is used."""
    use_plugin('simpleitk')
    yield
    reset_plugins()


def test_prefered_plugin():
    order = plugin_order()
    assert order["imread"][0] == "simpleitk"
    assert order["imsave"][0] == "simpleitk"
    assert order["imread_collection"][0] == "simpleitk"


def test_imread_as_gray():
    img = imread(testing.fetch('data/color.png'), as_gray=True)
    assert img.ndim == 2
    assert img.dtype == np.float64
    img = imread(testing.fetch('data/camera.png'), as_gray=True)
    # check that conversion does not happen for a gray image
    assert np.dtype(img.dtype).char in np.typecodes['AllInteger']


def test_bilevel():
    expected = np.zeros((10, 10))
    expected[::2] = 255

    img = imread(testing.fetch('data/checker_bilevel.png'))
    np.testing.assert_array_equal(img, expected)


def test_imread_truncated_jpg():
    with pytest.raises(RuntimeError):
        imread(testing.fetch('data/truncated.jpg'))


def test_imread_uint16():
    expected = np.load(testing.fetch('data/chessboard_GRAY_U8.npy'))
    img = imread(testing.fetch('data/chessboard_GRAY_U16.tif'))
    assert np.issubdtype(img.dtype, np.uint16)
    np.testing.assert_array_almost_equal(img, expected)


def test_imread_uint16_big_endian():
    expected = np.load(testing.fetch('data/chessboard_GRAY_U8.npy'))
    img = imread(testing.fetch('data/chessboard_GRAY_U16B.tif'), plugin="simpleitk")
    assert img.dtype.type == np.uint16
    np.testing.assert_array_almost_equal(img, expected)


@pytest.mark.parametrize("shape", [(10, 10), (10, 10, 3), (10, 10, 4)])
@pytest.mark.parametrize("dtype", [np.uint8, np.uint16, np.float32, np.float64])
# Neither NumPy's linspace nor called functions in our color module can
# really deal with the max and min value supported by float64. As we might
# deprecate this plugin system soon, just ignore it.
@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_imsave_roundtrip(shape, dtype, tmp_path):
    if np.issubdtype(dtype, np.floating):
        info_func = np.finfo
    else:
        info_func = np.iinfo
    expected = np.linspace(
        info_func(dtype).min,
        info_func(dtype).max,
        endpoint=True,
        num=np.prod(shape),
        dtype=dtype,
    )
    expected = expected.reshape(shape)
    file_path = tmp_path / "roundtrip.mha"
    imsave(file_path, expected)
    actual = imread(file_path)
    np.testing.assert_array_almost_equal(actual, expected)
