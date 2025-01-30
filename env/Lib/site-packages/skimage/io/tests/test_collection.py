import os
import itertools

import numpy as np
import imageio.v3 as iio3
from skimage import data_dir
from skimage.io.collection import ImageCollection, MultiImage, alphanumeric_key
from skimage.io import reset_plugins

from skimage._shared import testing
from skimage._shared.testing import assert_equal, assert_allclose, fetch

import pytest


try:
    has_pooch = True
except ModuleNotFoundError:
    has_pooch = False


def test_string_split():
    test_string = 'z23a'
    test_str_result = ['z', 23, 'a']
    assert_equal(alphanumeric_key(test_string), test_str_result)


def test_string_sort():
    filenames = [
        'f9.10.png',
        'f9.9.png',
        'f10.10.png',
        'f10.9.png',
        'e9.png',
        'e10.png',
        'em.png',
    ]
    expected_filenames = [
        'e9.png',
        'e10.png',
        'em.png',
        'f9.9.png',
        'f9.10.png',
        'f10.9.png',
        'f10.10.png',
    ]
    sorted_filenames = sorted(filenames, key=alphanumeric_key)
    assert_equal(expected_filenames, sorted_filenames)


def test_imagecollection_input():
    """Test function for ImageCollection. The new behavior (implemented
    in 0.16) allows the `pattern` argument to accept a list of strings
    as the input.

    Notes
    -----
        If correct, `images` will receive three images.
    """
    # Ensure that these images are part of the legacy datasets
    # this means they will always be available in the user's install
    # regardless of the availability of pooch
    pics = [
        fetch('data/coffee.png'),
        fetch('data/chessboard_GRAY.png'),
        fetch('data/rocket.jpg'),
    ]
    pattern = [os.path.join(data_dir, pic) for pic in pics]
    images = ImageCollection(pattern)
    assert len(images) == 3


class TestImageCollection:
    pics = [fetch('data/brick.png'), fetch('data/color.png'), fetch('data/moon.png')]
    pattern = pics[:2]
    pattern_same_shape = pics[::2]

    def setup_method(self):
        reset_plugins()
        # Generic image collection with images of different shapes.
        self.images = ImageCollection(self.pattern)
        # Image collection with images having shapes that match.
        self.images_matched = ImageCollection(self.pattern_same_shape)
        # Same images as a collection of frames
        self.frames_matched = MultiImage(self.pattern_same_shape)

    def test_len(self):
        assert len(self.images) == 2

    def test_getitem(self):
        num = len(self.images)
        for i in range(-num, num):
            assert isinstance(self.images[i], np.ndarray)
        assert_allclose(self.images[0], self.images[-num])

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

    @pytest.mark.skipif(not has_pooch, reason="needs pooch to download data")
    def test_custom_load_func_sequence(self):
        filename = fetch('data/no_time_for_that_tiny.gif')

        def reader(index):
            return iio3.imread(filename, index=index)

        ic = ImageCollection(range(24), load_func=reader)
        # the length of ic should be that of the given load_pattern sequence
        assert len(ic) == 24
        # GIF file has frames of size 25x14 with 4 channels (RGBA)
        assert ic[0].shape == (25, 14, 3)

    @pytest.mark.skipif(not has_pooch, reason="needs pooch to download data")
    def test_custom_load_func_w_kwarg(self):
        load_pattern = fetch('data/no_time_for_that_tiny.gif')

        def load_fn(f, step):
            vid = iio3.imiter(f)
            return list(itertools.islice(vid, None, None, step))

        ic = ImageCollection(load_pattern, load_func=load_fn, step=3)
        # Each file should map to one image (array).
        assert len(ic) == 1
        # GIF file has 24 frames, so 24 / 3 equals 8.
        assert len(ic[0]) == 8

    def test_custom_load_func(self):
        def load_fn(x):
            return x

        ic = ImageCollection(os.pathsep.join(self.pattern), load_func=load_fn)
        assert_equal(ic[0], self.pattern[0])

    def test_concatenate(self):
        array = self.images_matched.concatenate()
        expected_shape = (len(self.images_matched),) + self.images[0].shape
        assert_equal(array.shape, expected_shape)

    def test_concatenate_mismatched_image_shapes(self):
        with testing.raises(ValueError):
            self.images.concatenate()

    def test_multiimage_imagecollection(self):
        assert_equal(self.images_matched[0], self.frames_matched[0])
        assert_equal(self.images_matched[1], self.frames_matched[1])
