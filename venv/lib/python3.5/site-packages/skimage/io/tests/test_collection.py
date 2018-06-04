import os

import numpy as np
from skimage import data_dir
from skimage.io.collection import ImageCollection, alphanumeric_key
from skimage.io import reset_plugins

from skimage._shared import testing
from skimage._shared.testing import assert_equal, assert_allclose, TestCase


def test_string_split():
    test_string = 'z23a'
    test_str_result = ['z', 23, 'a']
    assert_equal(alphanumeric_key(test_string), test_str_result)


def test_string_sort():
    filenames = ['f9.10.png', 'f9.9.png', 'f10.10.png', 'f10.9.png',
                 'e9.png', 'e10.png', 'em.png']
    sorted_filenames = ['e9.png', 'e10.png', 'em.png', 'f9.9.png',
                        'f9.10.png', 'f10.9.png', 'f10.10.png']
    sorted_filenames = sorted(filenames, key=alphanumeric_key)
    assert_equal(sorted_filenames, sorted_filenames)


class TestImageCollection(TestCase):

    pattern = [os.path.join(data_dir, pic)
               for pic in ['camera.png', 'color.png']]

    pattern_matched = [os.path.join(data_dir, pic)
                       for pic in ['camera.png', 'moon.png']]

    @testing.fixture(autouse=True)
    def setUp(self):
        reset_plugins()
        # Generic image collection with images of different shapes.
        self.images = ImageCollection(self.pattern)
        # Image collection with images having shapes that match.
        self.images_matched = ImageCollection(self.pattern_matched)

    def test_len(self):
        assert len(self.images) == 2

    def test_getitem(self):
        num = len(self.images)
        for i in range(-num, num):
            assert type(self.images[i]) is np.ndarray
        assert_allclose(self.images[0],
                        self.images[-num])

        def return_img(n):
            return self.images[n]
        with testing.raises(IndexError):
            return_img(num)
        with testing.raises(IndexError):
            return_img(-num - 1)

    def test_slicing(self):
        assert type(self.images[:]) is ImageCollection
        assert len(self.images[:]) == 2
        assert len(self.images[:1]) == 1
        assert len(self.images[1:]) == 1
        assert_allclose(self.images[0], self.images[:1][0])
        assert_allclose(self.images[1], self.images[1:][0])
        assert_allclose(self.images[1], self.images[::-1][0])
        assert_allclose(self.images[0], self.images[::-1][1])

    def test_files_property(self):
        assert isinstance(self.images.files, list)

        def set_files(f):
            self.images.files = f
        with testing.raises(AttributeError):
            set_files('newfiles')

    def test_custom_load(self):
        load_pattern = [(1, 'one'), (2, 'two')]

        def load_fn(x):
            return x

        ic = ImageCollection(load_pattern, load_func=load_fn)
        assert_equal(ic[1], (2, 'two'))

    def test_custom_load_func(self):

        def load_fn(x):
            return x

        ic = ImageCollection(os.pathsep.join(self.pattern), load_func=load_fn)
        assert_equal(ic[0], self.pattern[0])

    def test_concatenate(self):
        array = self.images_matched.concatenate()
        expected_shape = (len(self.images_matched),) + self.images[0].shape
        assert_equal(array.shape, expected_shape)

    def test_concatentate_mismatched_image_shapes(self):
        with testing.raises(ValueError):
            self.images.concatenate()
