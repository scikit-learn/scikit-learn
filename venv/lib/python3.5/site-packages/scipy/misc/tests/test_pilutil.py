from __future__ import division, print_function, absolute_import

import os.path
import tempfile
import shutil
import numpy as np
import glob

import pytest
from pytest import raises as assert_raises
from numpy.testing import (assert_equal, assert_allclose,
                           assert_array_equal, assert_)
from scipy._lib._numpy_compat import suppress_warnings
from scipy import misc
from numpy.ma.testutils import assert_mask_equal

try:
    import PIL.Image
except ImportError:
    _have_PIL = False
else:
    _have_PIL = True


# Function / method decorator for skipping PIL tests on import failure
_pilskip = pytest.mark.skipif(not _have_PIL, reason='Need to import PIL for this test')

datapath = os.path.dirname(__file__)

@_pilskip
class TestPILUtil(object):
    def test_imresize(self):
        im = np.random.random((10, 20))
        for T in np.sctypes['float'] + [float]:
            # 1.1 rounds to below 1.1 for float16, 1.101 works
            with suppress_warnings() as sup:
                sup.filter(DeprecationWarning)
                im1 = misc.imresize(im, T(1.101))
            assert_equal(im1.shape, (11, 22))

    def test_imresize2(self):
        im = np.random.random((20, 30))
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im2 = misc.imresize(im, (30, 40), interp='bicubic')
        assert_equal(im2.shape, (30, 40))

    def test_imresize3(self):
        im = np.random.random((15, 30))
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im2 = misc.imresize(im, (30, 60), interp='nearest')
        assert_equal(im2.shape, (30, 60))

    def test_imresize4(self):
        im = np.array([[1, 2],
                       [3, 4]])
        # Check that resizing by target size, float and int are the same
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im2 = misc.imresize(im, (4, 4), mode='F')  # output size
            im3 = misc.imresize(im, 2., mode='F')  # fraction
            im4 = misc.imresize(im, 200, mode='F')  # percentage
        assert_equal(im2, im3)
        assert_equal(im2, im4)

    def test_imresize5(self):
        im = np.random.random((25, 15))
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im2 = misc.imresize(im, (30, 60), interp='lanczos')
        assert_equal(im2.shape, (30, 60))

    def test_bytescale(self):
        x = np.array([0, 1, 2], np.uint8)
        y = np.array([0, 1, 2])
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            assert_equal(misc.bytescale(x), x)
            assert_equal(misc.bytescale(y), [0, 128, 255])

    def test_bytescale_keywords(self):
        x = np.array([40, 60, 120, 200, 300, 500])
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            res_lowhigh = misc.bytescale(x, low=10, high=143)
            assert_equal(res_lowhigh, [10, 16, 33, 56, 85, 143])
            res_cmincmax = misc.bytescale(x, cmin=60, cmax=300)
            assert_equal(res_cmincmax, [0, 0, 64, 149, 255, 255])
            assert_equal(misc.bytescale(np.array([3, 3, 3]), low=4), [4, 4, 4])

    def test_bytescale_cscale_lowhigh(self):
        a = np.arange(10)
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            actual = misc.bytescale(a, cmin=3, cmax=6, low=100, high=200)
        expected = [100, 100, 100, 100, 133, 167, 200, 200, 200, 200]
        assert_equal(actual, expected)

    def test_bytescale_mask(self):
        a = np.ma.MaskedArray(data=[1, 2, 3], mask=[False, False, True])
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            actual = misc.bytescale(a)
        expected = [0, 255, 3]
        assert_equal(expected, actual)
        assert_mask_equal(a.mask, actual.mask)
        assert_(isinstance(actual, np.ma.MaskedArray))

    def test_bytescale_rounding(self):
        a = np.array([-0.5, 0.5, 1.5, 2.5, 3.5])
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            actual = misc.bytescale(a, cmin=0, cmax=10, low=0, high=10)
        expected = [0, 1, 2, 3, 4]
        assert_equal(actual, expected)

    def test_bytescale_low_greaterthan_high(self):
        with assert_raises(ValueError):
            with suppress_warnings() as sup:
                sup.filter(DeprecationWarning)
                misc.bytescale(np.arange(3), low=10, high=5)

    def test_bytescale_low_lessthan_0(self):
        with assert_raises(ValueError):
            with suppress_warnings() as sup:
                sup.filter(DeprecationWarning)
                misc.bytescale(np.arange(3), low=-1)

    def test_bytescale_high_greaterthan_255(self):
        with assert_raises(ValueError):
            with suppress_warnings() as sup:
                sup.filter(DeprecationWarning)
                misc.bytescale(np.arange(3), high=256)

    def test_bytescale_low_equals_high(self):
        a = np.arange(3)
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            actual = misc.bytescale(a, low=10, high=10)
        expected = [10, 10, 10]
        assert_equal(actual, expected)

    def test_imsave(self):
        picdir = os.path.join(datapath, "data")
        for png in glob.iglob(picdir + "/*.png"):
            with suppress_warnings() as sup:
                # PIL causes a Py3k ResourceWarning
                sup.filter(message="unclosed file")
                sup.filter(DeprecationWarning)
                img = misc.imread(png)
            tmpdir = tempfile.mkdtemp()
            try:
                fn1 = os.path.join(tmpdir, 'test.png')
                fn2 = os.path.join(tmpdir, 'testimg')
                with suppress_warnings() as sup:
                    # PIL causes a Py3k ResourceWarning
                    sup.filter(message="unclosed file")
                    sup.filter(DeprecationWarning)
                    misc.imsave(fn1, img)
                    misc.imsave(fn2, img, 'PNG')

                with suppress_warnings() as sup:
                    # PIL causes a Py3k ResourceWarning
                    sup.filter(message="unclosed file")
                    sup.filter(DeprecationWarning)
                    data1 = misc.imread(fn1)
                    data2 = misc.imread(fn2)
                assert_allclose(data1, img)
                assert_allclose(data2, img)
                assert_equal(data1.shape, img.shape)
                assert_equal(data2.shape, img.shape)
            finally:
                shutil.rmtree(tmpdir)


def check_fromimage(filename, irange, shape):
    fp = open(filename, "rb")
    with suppress_warnings() as sup:
        sup.filter(DeprecationWarning)
        img = misc.fromimage(PIL.Image.open(fp))
    fp.close()
    imin, imax = irange
    assert_equal(img.min(), imin)
    assert_equal(img.max(), imax)
    assert_equal(img.shape, shape)


@_pilskip
def test_fromimage():
    # Test generator for parametric tests
    # Tuples in the list are (filename, (datamin, datamax), shape).
    files = [('icon.png', (0, 255), (48, 48, 4)),
             ('icon_mono.png', (0, 255), (48, 48, 4)),
             ('icon_mono_flat.png', (0, 255), (48, 48, 3))]
    for fn, irange, shape in files:
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            check_fromimage(os.path.join(datapath, 'data', fn), irange, shape)


@_pilskip
def test_imread_indexed_png():
    # The file `foo3x5x4indexed.png` was created with this array
    # (3x5 is (height)x(width)):
    data = np.array([[[127, 0, 255, 255],
                      [127, 0, 255, 255],
                      [127, 0, 255, 255],
                      [127, 0, 255, 255],
                      [127, 0, 255, 255]],
                     [[192, 192, 255, 0],
                      [192, 192, 255, 0],
                      [0, 0, 255, 0],
                      [0, 0, 255, 0],
                      [0, 0, 255, 0]],
                     [[0, 31, 255, 255],
                      [0, 31, 255, 255],
                      [0, 31, 255, 255],
                      [0, 31, 255, 255],
                      [0, 31, 255, 255]]], dtype=np.uint8)

    filename = os.path.join(datapath, 'data', 'foo3x5x4indexed.png')
    with open(filename, 'rb') as f:
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im = misc.imread(f)
    assert_array_equal(im, data)


@_pilskip
def test_imread_1bit():
    # box1.png is a 48x48 grayscale image with bit depth 1.
    # The border pixels are 1 and the rest are 0.
    filename = os.path.join(datapath, 'data', 'box1.png')
    with open(filename, 'rb') as f:
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im = misc.imread(f)
    assert_equal(im.dtype, np.uint8)
    expected = np.zeros((48, 48), dtype=np.uint8)
    # When scaled up from 1 bit to 8 bits, 1 becomes 255.
    expected[:, 0] = 255
    expected[:, -1] = 255
    expected[0, :] = 255
    expected[-1, :] = 255
    assert_equal(im, expected)


@_pilskip
def test_imread_2bit():
    # blocks2bit.png is a 12x12 grayscale image with bit depth 2.
    # The pattern is 4 square subblocks of size 6x6.  Upper left
    # is all 0, upper right is all 1, lower left is all 2, lower
    # right is all 3.
    # When scaled up to 8 bits, the values become [0, 85, 170, 255].
    filename = os.path.join(datapath, 'data', 'blocks2bit.png')
    with open(filename, 'rb') as f:
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im = misc.imread(f)
    assert_equal(im.dtype, np.uint8)
    expected = np.zeros((12, 12), dtype=np.uint8)
    expected[:6, 6:] = 85
    expected[6:, :6] = 170
    expected[6:, 6:] = 255
    assert_equal(im, expected)


@_pilskip
def test_imread_4bit():
    # pattern4bit.png is a 12(h) x 31(w) grayscale image with bit depth 4.
    # The value in row j and column i is maximum(j, i) % 16.
    # When scaled up to 8 bits, the values become [0, 17, 34, ..., 255].
    filename = os.path.join(datapath, 'data', 'pattern4bit.png')
    with open(filename, 'rb') as f:
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            im = misc.imread(f)
    assert_equal(im.dtype, np.uint8)
    j, i = np.meshgrid(np.arange(12), np.arange(31), indexing='ij')
    expected = 17*(np.maximum(j, i) % 16).astype(np.uint8)
    assert_equal(im, expected)

