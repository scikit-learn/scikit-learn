import os.path
import numpy as np
import unittest

from tempfile import NamedTemporaryFile

from skimage import data_dir
from skimage.io import imread, imsave, use_plugin, reset_plugins
from skimage._shared import testing

try:
    import SimpleITK as sitk
    use_plugin('simpleitk')
except ImportError:
    sitk_available = False
else:
    sitk_available = True

np.random.seed(0)


def teardown():
    reset_plugins()


def setup_module(self):
    """The effect of the `plugin.use` call may be overridden by later imports.
    Call `use_plugin` directly before the tests to ensure that sitk is used.

    """
    try:
        use_plugin('simpleitk')
    except ImportError:
        pass


@testing.skipif(not sitk_available, reason="simpletk not installed")
def test_imread_flatten():
    # a color image is flattened
    img = imread(os.path.join(data_dir, 'color.png'), flatten=True)
    assert img.ndim == 2
    assert img.dtype == np.float64
    img = imread(os.path.join(data_dir, 'camera.png'), flatten=True)
    # check that flattening does not occur for an image that is grey already.
    assert np.sctype2char(img.dtype) in np.typecodes['AllInteger']


@testing.skipif(not sitk_available, reason="simpletk not installed")
def test_bilevel():
    expected = np.zeros((10, 10))
    expected[::2] = 255

    img = imread(os.path.join(data_dir, 'checker_bilevel.png'))
    np.testing.assert_array_equal(img, expected)

"""
#TODO: This test causes a Segmentation fault
@testing.skipif(not sitk_available)
def test_imread_truncated_jpg():
    assert_raises((RuntimeError, ValueError),
                  imread,
                  os.path.join(data_dir, 'truncated.jpg'))
"""


@testing.skipif(not sitk_available, reason="simpletk not installed")
def test_imread_uint16():
    expected = np.load(os.path.join(data_dir, 'chessboard_GRAY_U8.npy'))
    img = imread(os.path.join(data_dir, 'chessboard_GRAY_U16.tif'))
    assert np.issubdtype(img.dtype, np.uint16)
    np.testing.assert_array_almost_equal(img, expected)


@testing.skipif(not sitk_available, reason="simpletk not installed")
def test_imread_uint16_big_endian():
    expected = np.load(os.path.join(data_dir, 'chessboard_GRAY_U8.npy'))
    img = imread(os.path.join(data_dir, 'chessboard_GRAY_U16B.tif'))
    np.testing.assert_array_almost_equal(img, expected)


class TestSave(unittest.TestCase):
    def roundtrip(self, dtype, x):
        f = NamedTemporaryFile(suffix='.mha')
        fname = f.name
        f.close()
        imsave(fname, x)
        y = imread(fname)

        np.testing.assert_array_almost_equal(x, y)

    @testing.skipif(not sitk_available, reason="simpletk not installed")
    def test_imsave_roundtrip(self):
        for shape in [(10, 10), (10, 10, 3), (10, 10, 4)]:
            for dtype in (np.uint8, np.uint16, np.float32, np.float64):
                x = np.ones(shape, dtype=dtype) * np.random.rand(*shape)

                if np.issubdtype(dtype, np.floating):
                    yield self.roundtrip, dtype, x
                else:
                    x = (x * 255).astype(dtype)
                    yield self.roundtrip, dtype, x
