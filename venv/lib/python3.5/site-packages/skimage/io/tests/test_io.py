import os

import numpy as np
from skimage import io, data_dir

from skimage._shared import testing
from skimage._shared.testing import assert_array_equal


def test_stack_basic():
    x = np.arange(12).reshape(3, 4)
    io.push(x)

    assert_array_equal(io.pop(), x)


def test_stack_non_array():
    with testing.raises(ValueError):
        io.push([[1, 2, 3]])


def test_imread_url():
    # tweak data path so that file URI works on both unix and windows.
    data_path = data_dir.lstrip(os.path.sep)
    data_path = data_path.replace(os.path.sep, '/')
    image_url = 'file:///{0}/camera.png'.format(data_path)
    image = io.imread(image_url)
    assert image.shape == (512, 512)
