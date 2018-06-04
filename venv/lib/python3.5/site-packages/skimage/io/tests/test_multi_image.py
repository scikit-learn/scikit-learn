import os
import six

import numpy as np
from skimage import data_dir
from skimage.io import use_plugin
from skimage.io.collection import MultiImage, ImageCollection

from skimage._shared import testing
from skimage._shared.testing import assert_equal, assert_allclose, TestCase


class TestMultiImage(TestCase):
    @testing.fixture(autouse=True)
    def setUp(self):
        # This multipage TIF file was created with imagemagick:
        # convert im1.tif im2.tif -adjoin multipage.tif
        use_plugin('pil')
        paths = [os.path.join(data_dir, 'multipage_rgb.tif'),
                 os.path.join(data_dir, 'no_time_for_that_tiny.gif')]
        self.imgs = [MultiImage(paths[0]),
                     MultiImage(paths[0], conserve_memory=False),
                     MultiImage(paths[1]),
                     MultiImage(paths[1], conserve_memory=False),
                     ImageCollection(paths[0]),
                     ImageCollection(paths[1], conserve_memory=False),
                     ImageCollection(os.pathsep.join(paths))]

    def test_shapes(self):
        img = self.imgs[-1]
        imgs = img[:]
        assert imgs[0].shape == imgs[1].shape
        assert imgs[0].shape == (10, 10, 3)

    def test_len(self):
        assert len(self.imgs[0]) == len(self.imgs[1]) == 2
        assert len(self.imgs[2]) == len(self.imgs[3]) == 24
        assert len(self.imgs[4]) == 2
        assert len(self.imgs[5]) == 24
        assert len(self.imgs[6]) == 26, len(self.imgs[6])

    def test_slicing(self):
        img = self.imgs[-1]
        assert type(img[:]) is ImageCollection
        assert len(img[:]) == 26, len(img[:])
        assert len(img[:1]) == 1
        assert len(img[1:]) == 25
        assert_allclose(img[0], img[:1][0])
        assert_allclose(img[1], img[1:][0])
        assert_allclose(img[-1], img[::-1][0])
        assert_allclose(img[0], img[::-1][-1])

    def test_getitem(self):
        for img in self.imgs:
            num = len(img)

            for i in range(-num, num):
                assert type(img[i]) is np.ndarray
            assert_allclose(img[0], img[-num])

            with testing.raises(AssertionError):
                assert_allclose(img[0], img[1])

            with testing.raises(IndexError):
                img[num]
            with testing.raises(IndexError):
                img[-num - 1]

    def test_files_property(self):
        for img in self.imgs:
            if isinstance(img, ImageCollection):
                continue

            assert isinstance(img.filename, six.string_types)

            with testing.raises(AttributeError):
                img.filename = "newfile"

    def test_conserve_memory_property(self):
        for img in self.imgs:
            assert isinstance(img.conserve_memory, bool)

            with testing.raises(AttributeError):
                img.conserve_memory = True

    def test_concatenate(self):
        for img in self.imgs:
            if img[0].shape != img[-1].shape:
                with testing.raises(ValueError):
                    img.concatenate()
                continue
            array = img.concatenate()
            assert_equal(array.shape, (len(img),) + img[0].shape)
