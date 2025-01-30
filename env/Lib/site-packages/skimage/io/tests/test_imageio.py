from tempfile import NamedTemporaryFile

import numpy as np
from skimage.io import imread, imsave, plugin_order

from skimage._shared import testing
from skimage._shared.testing import fetch, assert_stacklevel

import pytest


def test_prefered_plugin():
    # Don't call use_plugin("imageio") before, this way we test that imageio is used
    # by default
    order = plugin_order()
    assert order["imread"][0] == "imageio"
    assert order["imsave"][0] == "imageio"
    assert order["imread_collection"][0] == "imageio"


def test_imageio_as_gray():
    img = imread(fetch('data/color.png'), as_gray=True)
    assert img.ndim == 2
    assert img.dtype == np.float64
    img = imread(fetch('data/camera.png'), as_gray=True)
    # check that conversion does not happen for a gray image
    assert np.dtype(img.dtype).char in np.typecodes['AllInteger']


def test_imageio_palette():
    img = imread(fetch('data/palette_color.png'))
    assert img.ndim == 3


def test_imageio_truncated_jpg():
    # imageio>2.0 uses Pillow / PIL to try and load the file.
    # Oddly, PIL explicitly raises a SyntaxError when the file read fails.
    # The exception type changed from SyntaxError to OSError in PIL 8.2.0, so
    # allow for either to be raised.
    with testing.raises((OSError, SyntaxError)):
        imread(fetch('data/truncated.jpg'))


class TestSave:
    @pytest.mark.parametrize(
        "shape,dtype",
        [
            # float32, float64 can't be saved as PNG and raise
            # uint32 is not roundtripping properly
            ((10, 10), np.uint8),
            ((10, 10), np.uint16),
            ((10, 10, 2), np.uint8),
            ((10, 10, 3), np.uint8),
            ((10, 10, 4), np.uint8),
        ],
    )
    def test_imsave_roundtrip(self, shape, dtype, tmp_path):
        if np.issubdtype(dtype, np.floating):
            min_ = 0
            max_ = 1
        else:
            min_ = 0
            max_ = np.iinfo(dtype).max
        expected = np.linspace(
            min_, max_, endpoint=True, num=np.prod(shape), dtype=dtype
        )
        expected = expected.reshape(shape)
        file_path = tmp_path / "roundtrip.png"
        imsave(file_path, expected)
        actual = imread(file_path)
        np.testing.assert_array_almost_equal(actual, expected)

    def test_bool_array_save(self):
        with NamedTemporaryFile(suffix='.png') as f:
            fname = f.name

        with pytest.warns(UserWarning, match=r'.* is a boolean image') as record:
            a = np.zeros((5, 5), bool)
            a[2, 2] = True
            imsave(fname, a)
        assert_stacklevel(record)


def test_return_class():
    testing.assert_equal(type(imread(fetch('data/color.png'))), np.ndarray)
