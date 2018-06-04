import os
from tempfile import NamedTemporaryFile

import numpy as np
from skimage import data_dir
from skimage.io import imread, imsave, use_plugin, reset_plugins

from skimage._shared import testing
from skimage._shared.testing import assert_array_almost_equal, TestCase


try:
    import imageio as _imageio
except ImportError:
    imageio_available = False
else:
    imageio_available = True


def setup():
    if imageio_available:
        np.random.seed(0)
        use_plugin('imageio')


def teardown():
    reset_plugins()


@testing.skipif(not imageio_available, reason="imageio not installed")
def test_imageio_flatten():
    # a color image is flattened
    img = imread(os.path.join(data_dir, 'color.png'), flatten=True)
    assert img.ndim == 2
    assert img.dtype == np.float64
    img = imread(os.path.join(data_dir, 'camera.png'), flatten=True)
    # check that flattening does not occur for an image that is grey already.
    assert np.sctype2char(img.dtype) in np.typecodes['AllInteger']


@testing.skipif(not imageio_available, reason="imageio not installed")
def test_imageio_palette():
    img = imread(os.path.join(data_dir, 'palette_color.png'))
    assert img.ndim == 3


@testing.skipif(not imageio_available, reason="imageio not installed")
def test_imageio_truncated_jpg():
    # imageio>2.0 uses Pillow / PIL to try and load the file.
    # Oddly, PIL explicitly raises a SyntaxError when the file read fails.
    with testing.raises(SyntaxError):
        imread(os.path.join(data_dir, 'truncated.jpg'))


class TestSave(TestCase):

    def roundtrip(self, x, scaling=1):
        f = NamedTemporaryFile(suffix='.png')
        fname = f.name
        f.close()
        imsave(fname, x)
        y = imread(fname)

        assert_array_almost_equal((x * scaling).astype(np.int32), y)

    @testing.skipif(not imageio_available, reason="imageio not installed")
    def test_imsave_roundtrip(self):
        dtype = np.uint8
        for shape in [(10, 10), (10, 10, 3), (10, 10, 4)]:
            x = np.ones(shape, dtype=dtype) * np.random.rand(*shape)

            if np.issubdtype(dtype, np.floating):
                yield self.roundtrip, x, 255
            else:
                x = (x * 255).astype(dtype)
                yield self.roundtrip, x
