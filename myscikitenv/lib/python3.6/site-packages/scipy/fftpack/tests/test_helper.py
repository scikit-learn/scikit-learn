# Created by Pearu Peterson, September 2002

from __future__ import division, print_function, absolute_import

__usage__ = """
Build fftpack:
  python setup_fftpack.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.fftpack.test(<level>)'
Run tests if fftpack is not installed:
  python tests/test_helper.py [<level>]
"""

from pytest import raises as assert_raises
from numpy.testing import assert_array_almost_equal, assert_equal, assert_
from scipy.fftpack import fftshift,ifftshift,fftfreq,rfftfreq
from scipy.fftpack.helper import (next_fast_len,
                                  _init_nd_shape_and_axes,
                                  _init_nd_shape_and_axes_sorted)

from numpy import pi, random
import numpy as np

class TestFFTShift(object):

    def test_definition(self):
        x = [0,1,2,3,4,-4,-3,-2,-1]
        y = [-4,-3,-2,-1,0,1,2,3,4]
        assert_array_almost_equal(fftshift(x),y)
        assert_array_almost_equal(ifftshift(y),x)
        x = [0,1,2,3,4,-5,-4,-3,-2,-1]
        y = [-5,-4,-3,-2,-1,0,1,2,3,4]
        assert_array_almost_equal(fftshift(x),y)
        assert_array_almost_equal(ifftshift(y),x)

    def test_inverse(self):
        for n in [1,4,9,100,211]:
            x = random.random((n,))
            assert_array_almost_equal(ifftshift(fftshift(x)),x)


class TestFFTFreq(object):

    def test_definition(self):
        x = [0,1,2,3,4,-4,-3,-2,-1]
        assert_array_almost_equal(9*fftfreq(9),x)
        assert_array_almost_equal(9*pi*fftfreq(9,pi),x)
        x = [0,1,2,3,4,-5,-4,-3,-2,-1]
        assert_array_almost_equal(10*fftfreq(10),x)
        assert_array_almost_equal(10*pi*fftfreq(10,pi),x)


class TestRFFTFreq(object):

    def test_definition(self):
        x = [0,1,1,2,2,3,3,4,4]
        assert_array_almost_equal(9*rfftfreq(9),x)
        assert_array_almost_equal(9*pi*rfftfreq(9,pi),x)
        x = [0,1,1,2,2,3,3,4,4,5]
        assert_array_almost_equal(10*rfftfreq(10),x)
        assert_array_almost_equal(10*pi*rfftfreq(10,pi),x)


class TestNextOptLen(object):

    def test_next_opt_len(self):
        random.seed(1234)

        def nums():
            for j in range(1, 1000):
                yield j
            yield 2**5 * 3**5 * 4**5 + 1

        for n in nums():
            m = next_fast_len(n)
            msg = "n=%d, m=%d" % (n, m)

            assert_(m >= n, msg)

            # check regularity
            k = m
            for d in [2, 3, 5]:
                while True:
                    a, b = divmod(k, d)
                    if b == 0:
                        k = a
                    else:
                        break
            assert_equal(k, 1, err_msg=msg)

    def test_np_integers(self):
        ITYPES = [np.int16, np.int32, np.int64, np.uint16, np.uint32, np.uint64]
        for ityp in ITYPES:
            x = ityp(12345)
            testN = next_fast_len(x)
            assert_equal(testN, next_fast_len(int(x)))

    def test_next_opt_len_strict(self):
        hams = {
            1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 8, 8: 8, 14: 15, 15: 15,
            16: 16, 17: 18, 1021: 1024, 1536: 1536, 51200000: 51200000,
            510183360: 510183360, 510183360 + 1: 512000000,
            511000000: 512000000,
            854296875: 854296875, 854296875 + 1: 859963392,
            196608000000: 196608000000, 196608000000 + 1: 196830000000,
            8789062500000: 8789062500000, 8789062500000 + 1: 8796093022208,
            206391214080000: 206391214080000,
            206391214080000 + 1: 206624260800000,
            470184984576000: 470184984576000,
            470184984576000 + 1: 470715894135000,
            7222041363087360: 7222041363087360,
            7222041363087360 + 1: 7230196133913600,
            # power of 5    5**23
            11920928955078125: 11920928955078125,
            11920928955078125 - 1: 11920928955078125,
            # power of 3    3**34
            16677181699666569: 16677181699666569,
            16677181699666569 - 1: 16677181699666569,
            # power of 2   2**54
            18014398509481984: 18014398509481984,
            18014398509481984 - 1: 18014398509481984,
            # above this, int(ceil(n)) == int(ceil(n+1))
            19200000000000000: 19200000000000000,
            19200000000000000 + 1: 19221679687500000,
            288230376151711744: 288230376151711744,
            288230376151711744 + 1: 288325195312500000,
            288325195312500000 - 1: 288325195312500000,
            288325195312500000: 288325195312500000,
            288325195312500000 + 1: 288555831593533440,
            # power of 3    3**83
            3990838394187339929534246675572349035227 - 1:
                3990838394187339929534246675572349035227,
            3990838394187339929534246675572349035227:
                3990838394187339929534246675572349035227,
            # power of 2     2**135
            43556142965880123323311949751266331066368 - 1:
                43556142965880123323311949751266331066368,
            43556142965880123323311949751266331066368:
                43556142965880123323311949751266331066368,
            # power of 5      5**57
            6938893903907228377647697925567626953125 - 1:
                6938893903907228377647697925567626953125,
            6938893903907228377647697925567626953125:
                6938893903907228377647697925567626953125,
            # http://www.drdobbs.com/228700538
            # 2**96 * 3**1 * 5**13
            290142196707511001929482240000000000000 - 1:
                290142196707511001929482240000000000000,
            290142196707511001929482240000000000000:
                290142196707511001929482240000000000000,
            290142196707511001929482240000000000000 + 1:
                290237644800000000000000000000000000000,
            # 2**36 * 3**69 * 5**7
            4479571262811807241115438439905203543080960000000 - 1:
                4479571262811807241115438439905203543080960000000,
            4479571262811807241115438439905203543080960000000:
                4479571262811807241115438439905203543080960000000,
            4479571262811807241115438439905203543080960000000 + 1:
                4480327901140333639941336854183943340032000000000,
            # 2**37 * 3**44 * 5**42
            30774090693237851027531250000000000000000000000000000000000000 - 1:
                30774090693237851027531250000000000000000000000000000000000000,
            30774090693237851027531250000000000000000000000000000000000000:
                30774090693237851027531250000000000000000000000000000000000000,
            30774090693237851027531250000000000000000000000000000000000000 + 1:
                30778180617309082445871527002041377406962596539492679680000000,
        }
        for x, y in hams.items():
            assert_equal(next_fast_len(x), y)


class Test_init_nd_shape_and_axes(object):

    def test_py_0d_defaults(self):
        x = 4
        shape = None
        axes = None

        shape_expected = np.array([])
        axes_expected = np.array([])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_0d_defaults(self):
        x = np.array(7.)
        shape = None
        axes = None

        shape_expected = np.array([])
        axes_expected = np.array([])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_py_1d_defaults(self):
        x = [1, 2, 3]
        shape = None
        axes = None

        shape_expected = np.array([3])
        axes_expected = np.array([0])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_1d_defaults(self):
        x = np.arange(0, 1, .1)
        shape = None
        axes = None

        shape_expected = np.array([10])
        axes_expected = np.array([0])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_py_2d_defaults(self):
        x = [[1, 2, 3, 4],
             [5, 6, 7, 8]]
        shape = None
        axes = None

        shape_expected = np.array([2, 4])
        axes_expected = np.array([0, 1])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_2d_defaults(self):
        x = np.arange(0, 1, .1).reshape(5, 2)
        shape = None
        axes = None

        shape_expected = np.array([5, 2])
        axes_expected = np.array([0, 1])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_5d_defaults(self):
        x = np.zeros([6, 2, 5, 3, 4])
        shape = None
        axes = None

        shape_expected = np.array([6, 2, 5, 3, 4])
        axes_expected = np.array([0, 1, 2, 3, 4])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_5d_set_shape(self):
        x = np.zeros([6, 2, 5, 3, 4])
        shape = [10, -1, -1, 1, 4]
        axes = None

        shape_expected = np.array([10, 2, 5, 1, 4])
        axes_expected = np.array([0, 1, 2, 3, 4])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_5d_set_axes(self):
        x = np.zeros([6, 2, 5, 3, 4])
        shape = None
        axes = [4, 1, 2]

        shape_expected = np.array([4, 2, 5])
        axes_expected = np.array([4, 1, 2])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_5d_set_axes_sorted(self):
        x = np.zeros([6, 2, 5, 3, 4])
        shape = None
        axes = [4, 1, 2]

        shape_expected = np.array([2, 5, 4])
        axes_expected = np.array([1, 2, 4])

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_5d_set_shape_axes(self):
        x = np.zeros([6, 2, 5, 3, 4])
        shape = [10, -1, 2]
        axes = [1, 0, 3]

        shape_expected = np.array([10, 6, 2])
        axes_expected = np.array([1, 0, 3])

        shape_res, axes_res = _init_nd_shape_and_axes(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_np_5d_set_shape_axes_sorted(self):
        x = np.zeros([6, 2, 5, 3, 4])
        shape = [10, -1, 2]
        axes = [1, 0, 3]

        shape_expected = np.array([6, 10, 2])
        axes_expected = np.array([0, 1, 3])

        shape_res, axes_res = _init_nd_shape_and_axes_sorted(x, shape, axes)

        assert_equal(shape_res, shape_expected)
        assert_equal(axes_res, axes_expected)

    def test_errors(self):
        with assert_raises(ValueError,
                           match="when given, axes values must be a scalar"
                           " or vector"):
            _init_nd_shape_and_axes([0], shape=None, axes=[[1, 2], [3, 4]])

        with assert_raises(ValueError,
                           match="when given, axes values must be integers"):
            _init_nd_shape_and_axes([0], shape=None, axes=[1., 2., 3., 4.])

        with assert_raises(ValueError,
                           match="axes exceeds dimensionality of input"):
            _init_nd_shape_and_axes([0], shape=None, axes=[1])

        with assert_raises(ValueError,
                           match="axes exceeds dimensionality of input"):
            _init_nd_shape_and_axes([0], shape=None, axes=[-2])

        with assert_raises(ValueError,
                           match="all axes must be unique"):
            _init_nd_shape_and_axes([0], shape=None, axes=[0, 0])

        with assert_raises(ValueError,
                           match="when given, shape values must be a scalar "
                           "or vector"):
            _init_nd_shape_and_axes([0], shape=[[1, 2], [3, 4]], axes=None)

        with assert_raises(ValueError,
                           match="when given, shape values must be integers"):
            _init_nd_shape_and_axes([0], shape=[1., 2., 3., 4.], axes=None)

        with assert_raises(ValueError,
                           match="when given, axes and shape arguments"
                           " have to be of the same length"):
            _init_nd_shape_and_axes(np.zeros([1, 1, 1, 1]),
                                    shape=[1, 2, 3], axes=[1])

        with assert_raises(ValueError,
                           match="invalid number of data points"
                           r" \(\[0\]\) specified"):
            _init_nd_shape_and_axes([0], shape=[0], axes=None)

        with assert_raises(ValueError,
                           match="invalid number of data points"
                           r" \(\[-2\]\) specified"):
            _init_nd_shape_and_axes([0], shape=-2, axes=None)
