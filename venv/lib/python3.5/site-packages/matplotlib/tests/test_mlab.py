from __future__ import absolute_import, division, print_function

import six

import tempfile
import warnings

from numpy.testing import (assert_allclose, assert_almost_equal,
                           assert_array_equal, assert_array_almost_equal_nulp)
import numpy.ma.testutils as matest
import numpy as np
import datetime as datetime
import pytest

import matplotlib.mlab as mlab
import matplotlib.cbook as cbook
from matplotlib.cbook.deprecation import MatplotlibDeprecationWarning


try:
    from mpl_toolkits.natgrid import _natgrid
    HAS_NATGRID = True
except ImportError:
    HAS_NATGRID = False


'''
A lot of mlab.py has been deprecated in Matplotlib 2.2 and is scheduled for
removal in the future. The tests that use deprecated methods have a block
to catch the deprecation warning, and can be removed with the mlab code is
removed.
'''


def test_colinear_pca():
    with pytest.warns(MatplotlibDeprecationWarning):
        a = mlab.PCA._get_colinear()
        pca = mlab.PCA(a)

    assert_allclose(pca.fracs[2:], 0., atol=1e-8)
    assert_allclose(pca.Y[:, 2:], 0., atol=1e-8)


@pytest.mark.parametrize('input', [
    # test odd lengths
    [1, 2, 3],
    # test even lengths
    [1, 2, 3, 4],
    # derived from email sent by jason-sage to MPL-user on 20090914
    [1, 1, 2, 2, 1, 2, 4, 3, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 7, 6, 4, 5, 5],
],
ids=[
    'odd length',
    'even length',
    'custom data',
])
@pytest.mark.parametrize('percentile', [
    0,
    50,
    75,
    100,
    [0, 75, 100],
])
def test_prctile(input, percentile):
    with pytest.warns(MatplotlibDeprecationWarning):
        assert_allclose(mlab.prctile(input, percentile),
                        np.percentile(input, percentile))


@pytest.mark.parametrize('xmin, xmax, N', [
    (.01, 1000., 6),
    (.03, 1313., 7),
    (.03, 1313., 0),
    (.03, 1313., 1),
], ids=[
    'tens',
    'primes',
    'none',
    'single',
])
def test_logspace(xmin, xmax, N):
    with pytest.warns(MatplotlibDeprecationWarning):
        res = mlab.logspace(xmin, xmax, N)
    targ = np.logspace(np.log10(xmin), np.log10(xmax), N)
    assert_allclose(targ, res)
    assert res.size == N


class TestStride(object):
    def get_base(self, x):
        y = x
        while y.base is not None:
            y = y.base
        return y

    def calc_window_target(self, x, NFFT, noverlap=0, axis=0):
        '''This is an adaptation of the original window extraction
        algorithm.  This is here to test to make sure the new implementation
        has the same result'''
        step = NFFT - noverlap
        ind = np.arange(0, len(x) - NFFT + 1, step)
        n = len(ind)
        result = np.zeros((NFFT, n))

        # do the ffts of the slices
        for i in range(n):
            result[:, i] = x[ind[i]:ind[i]+NFFT]
        if axis == 1:
            result = result.T
        return result

    @pytest.mark.parametrize('shape', [(), (10, 1)], ids=['0D', '2D'])
    def test_stride_windows_invalid_input_shape(self, shape):
        x = np.arange(np.prod(shape)).reshape(shape)
        with pytest.raises(ValueError):
            mlab.stride_windows(x, 5)

    @pytest.mark.parametrize('n, noverlap',
                             [(0, None), (11, None), (2, 2), (2, 3)],
                             ids=['n less than 1', 'n greater than input',
                                  'noverlap greater than n',
                                  'noverlap equal to n'])
    def test_stride_windows_invalid_params(self, n, noverlap):
        x = np.arange(10)
        with pytest.raises(ValueError):
            mlab.stride_windows(x, n, noverlap)

    @pytest.mark.parametrize('shape', [(), (10, 1)], ids=['0D', '2D'])
    def test_stride_repeat_invalid_input_shape(self, shape):
        x = np.arange(np.prod(shape)).reshape(shape)
        with pytest.raises(ValueError):
            mlab.stride_repeat(x, 5)

    @pytest.mark.parametrize('axis', [-1, 2],
                             ids=['axis less than 0',
                                  'axis greater than input shape'])
    def test_stride_repeat_invalid_axis(self, axis):
        x = np.array(0)
        with pytest.raises(ValueError):
            mlab.stride_repeat(x, 5, axis=axis)

    def test_stride_repeat_n_lt_1_ValueError(self):
        x = np.arange(10)
        with pytest.raises(ValueError):
            mlab.stride_repeat(x, 0)

    @pytest.mark.parametrize('axis', [0, 1], ids=['axis0', 'axis1'])
    @pytest.mark.parametrize('n', [1, 5], ids=['n1', 'n5'])
    def test_stride_repeat(self, n, axis):
        x = np.arange(10)
        y = mlab.stride_repeat(x, n, axis=axis)

        expected_shape = [10, 10]
        expected_shape[axis] = n
        yr = np.repeat(np.expand_dims(x, axis), n, axis=axis)

        assert yr.shape == y.shape
        assert_array_equal(yr, y)
        assert tuple(expected_shape) == y.shape
        assert self.get_base(y) is x

    @pytest.mark.parametrize('axis', [0, 1], ids=['axis0', 'axis1'])
    @pytest.mark.parametrize('n, noverlap',
                             [(1, 0), (5, 0), (15, 2), (13, -3)],
                             ids=['n1-noverlap0', 'n5-noverlap0',
                                  'n15-noverlap2', 'n13-noverlapn3'])
    def test_stride_windows(self, n, noverlap, axis):
        x = np.arange(100)
        y = mlab.stride_windows(x, n, noverlap=noverlap, axis=axis)

        expected_shape = [0, 0]
        expected_shape[axis] = n
        expected_shape[1 - axis] = 100 // (n - noverlap)
        yt = self.calc_window_target(x, n, noverlap=noverlap, axis=axis)

        assert yt.shape == y.shape
        assert_array_equal(yt, y)
        assert tuple(expected_shape) == y.shape
        assert self.get_base(y) is x

    @pytest.mark.parametrize('axis', [0, 1], ids=['axis0', 'axis1'])
    def test_stride_windows_n32_noverlap0_unflatten(self, axis):
        n = 32
        x = np.arange(n)[np.newaxis]
        x1 = np.tile(x, (21, 1))
        x2 = x1.flatten()
        y = mlab.stride_windows(x2, n, axis=axis)

        if axis == 0:
            x1 = x1.T
        assert y.shape == x1.shape
        assert_array_equal(y, x1)

    def test_stride_ensure_integer_type(self):
        N = 100
        x = np.empty(N + 20, dtype='>f4')
        x.fill(np.NaN)
        y = x[10:-10]
        y.fill(0.3)
        # previous to #3845 lead to corrupt access
        y_strided = mlab.stride_windows(y, n=33, noverlap=0.6)
        assert_array_equal(y_strided, 0.3)
        # previous to #3845 lead to corrupt access
        y_strided = mlab.stride_windows(y, n=33.3, noverlap=0)
        assert_array_equal(y_strided, 0.3)
        # even previous to #3845 could not find any problematic
        # configuration however, let's be sure it's not accidentally
        # introduced
        y_strided = mlab.stride_repeat(y, n=33.815)
        assert_array_equal(y_strided, 0.3)


@pytest.fixture
def tempcsv():
    if six.PY2:
        fd = tempfile.TemporaryFile(suffix='csv', mode="wb+")
    else:
        fd = tempfile.TemporaryFile(suffix='csv', mode="w+", newline='')
    with fd:
        yield fd


def test_recarray_csv_roundtrip(tempcsv):
    expected = np.recarray((99,),
                           [(str('x'), float),
                            (str('y'), float),
                            (str('t'), float)])
    # initialising all values: uninitialised memory sometimes produces
    # floats that do not round-trip to string and back.
    expected['x'][:] = np.linspace(-1e9, -1, 99)
    expected['y'][:] = np.linspace(1, 1e9, 99)
    expected['t'][:] = np.linspace(0, 0.01, 99)
    with pytest.warns(MatplotlibDeprecationWarning):
        mlab.rec2csv(expected, tempcsv)
        tempcsv.seek(0)
        actual = mlab.csv2rec(tempcsv)

    assert_allclose(expected['x'], actual['x'])
    assert_allclose(expected['y'], actual['y'])
    assert_allclose(expected['t'], actual['t'])


def test_rec2csv_bad_shape_ValueError(tempcsv):
    bad = np.recarray((99, 4), [(str('x'), float),
                                (str('y'), float)])

    # the bad recarray should trigger a ValueError for having ndim > 1.
    with pytest.warns(MatplotlibDeprecationWarning):
        with pytest.raises(ValueError):
            mlab.rec2csv(bad, tempcsv)


def test_csv2rec_names_with_comments(tempcsv):

    tempcsv.write('# comment\n1,2,3\n4,5,6\n')
    tempcsv.seek(0)
    with pytest.warns(MatplotlibDeprecationWarning):
        array = mlab.csv2rec(tempcsv, names='a,b,c')
    assert len(array) == 2
    assert len(array.dtype) == 3


@pytest.mark.parametrize('input, kwargs', [
    ('01/11/14\n'
     '03/05/76 12:00:01 AM\n'
     '07/09/83 5:17:34 PM\n'
     '06/20/2054 2:31:45 PM\n'
     '10/31/00 11:50:23 AM\n',
     {}),
    ('11/01/14\n'
     '05/03/76 12:00:01 AM\n'
     '09/07/83 5:17:34 PM\n'
     '20/06/2054 2:31:45 PM\n'
     '31/10/00 11:50:23 AM\n',
     {'dayfirst': True}),
    ('14/01/11\n'
     '76/03/05 12:00:01 AM\n'
     '83/07/09 5:17:34 PM\n'
     '2054/06/20 2:31:45 PM\n'
     '00/10/31 11:50:23 AM\n',
     {'yearfirst': True}),
], ids=['usdate', 'dayfirst', 'yearfirst'])
def test_csv2rec_dates(tempcsv, input, kwargs):
    tempcsv.write(input)
    expected = [datetime.datetime(2014, 1, 11, 0, 0),
                datetime.datetime(1976, 3, 5, 0, 0, 1),
                datetime.datetime(1983, 7, 9, 17, 17, 34),
                datetime.datetime(2054, 6, 20, 14, 31, 45),
                datetime.datetime(2000, 10, 31, 11, 50, 23)]
    tempcsv.seek(0)
    with pytest.warns(MatplotlibDeprecationWarning):
        array = mlab.csv2rec(tempcsv, names='a', **kwargs)
    assert_array_equal(array['a'].tolist(), expected)


def test_rec2txt_basic():
    # str() calls around field names necessary b/c as of numpy 1.11
    # dtype doesn't like unicode names (caused by unicode_literals import)
    a = np.array([(1.0, 2, 'foo', 'bing'),
                  (2.0, 3, 'bar', 'blah')],
                 dtype=np.dtype([(str('x'), np.float32),
                                 (str('y'), np.int8),
                                 (str('s'), str, 3),
                                 (str('s2'), str, 4)]))
    truth = ('       x   y   s     s2\n'
             '   1.000   2   foo   bing   \n'
             '   2.000   3   bar   blah   ').splitlines()
    with pytest.warns(MatplotlibDeprecationWarning):
        assert mlab.rec2txt(a).splitlines() == truth


class TestWindow(object):
    def setup(self):
        np.random.seed(0)
        n = 1000

        self.sig_rand = np.random.standard_normal(n) + 100.
        self.sig_ones = np.ones(n)

    def check_window_apply_repeat(self, x, window, NFFT, noverlap):
        '''This is an adaptation of the original window application
        algorithm.  This is here to test to make sure the new implementation
        has the same result'''
        step = NFFT - noverlap
        ind = np.arange(0, len(x) - NFFT + 1, step)
        n = len(ind)
        result = np.zeros((NFFT, n))

        if cbook.iterable(window):
            windowVals = window
        else:
            windowVals = window(np.ones((NFFT,), x.dtype))

        # do the ffts of the slices
        for i in range(n):
            result[:, i] = windowVals * x[ind[i]:ind[i]+NFFT]
        return result

    def test_window_none_rand(self):
        res = mlab.window_none(self.sig_ones)
        assert_array_equal(res, self.sig_ones)

    def test_window_none_ones(self):
        res = mlab.window_none(self.sig_rand)
        assert_array_equal(res, self.sig_rand)

    def test_window_hanning_rand(self):
        targ = np.hanning(len(self.sig_rand)) * self.sig_rand
        res = mlab.window_hanning(self.sig_rand)

        assert_allclose(targ, res, atol=1e-06)

    def test_window_hanning_ones(self):
        targ = np.hanning(len(self.sig_ones))
        res = mlab.window_hanning(self.sig_ones)

        assert_allclose(targ, res, atol=1e-06)

    def test_apply_window_1D_axis1_ValueError(self):
        x = self.sig_rand
        window = mlab.window_hanning
        with pytest.raises(ValueError):
            mlab.apply_window(x, window, axis=1, return_window=False)

    def test_apply_window_1D_els_wrongsize_ValueError(self):
        x = self.sig_rand
        window = mlab.window_hanning(np.ones(x.shape[0]-1))
        with pytest.raises(ValueError):
            mlab.apply_window(x, window)

    def test_apply_window_0D_ValueError(self):
        x = np.array(0)
        window = mlab.window_hanning
        with pytest.raises(ValueError):
            mlab.apply_window(x, window, axis=1, return_window=False)

    def test_apply_window_3D_ValueError(self):
        x = self.sig_rand[np.newaxis][np.newaxis]
        window = mlab.window_hanning
        with pytest.raises(ValueError):
            mlab.apply_window(x, window, axis=1, return_window=False)

    def test_apply_window_hanning_1D(self):
        x = self.sig_rand
        window = mlab.window_hanning
        window1 = mlab.window_hanning(np.ones(x.shape[0]))
        y, window2 = mlab.apply_window(x, window, return_window=True)
        yt = window(x)
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)
        assert_array_equal(window1, window2)

    def test_apply_window_hanning_1D_axis0(self):
        x = self.sig_rand
        window = mlab.window_hanning
        y = mlab.apply_window(x, window, axis=0, return_window=False)
        yt = window(x)
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_hanning_els_1D_axis0(self):
        x = self.sig_rand
        window = mlab.window_hanning(np.ones(x.shape[0]))
        window1 = mlab.window_hanning
        y = mlab.apply_window(x, window, axis=0, return_window=False)
        yt = window1(x)
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_hanning_2D_axis0(self):
        x = np.random.standard_normal([1000, 10]) + 100.
        window = mlab.window_hanning
        y = mlab.apply_window(x, window, axis=0, return_window=False)
        yt = np.zeros_like(x)
        for i in range(x.shape[1]):
            yt[:, i] = window(x[:, i])
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_hanning_els1_2D_axis0(self):
        x = np.random.standard_normal([1000, 10]) + 100.
        window = mlab.window_hanning(np.ones(x.shape[0]))
        window1 = mlab.window_hanning
        y = mlab.apply_window(x, window, axis=0, return_window=False)
        yt = np.zeros_like(x)
        for i in range(x.shape[1]):
            yt[:, i] = window1(x[:, i])
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_hanning_els2_2D_axis0(self):
        x = np.random.standard_normal([1000, 10]) + 100.
        window = mlab.window_hanning
        window1 = mlab.window_hanning(np.ones(x.shape[0]))
        y, window2 = mlab.apply_window(x, window, axis=0, return_window=True)
        yt = np.zeros_like(x)
        for i in range(x.shape[1]):
            yt[:, i] = window1*x[:, i]
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)
        assert_array_equal(window1, window2)

    def test_apply_window_hanning_els3_2D_axis0(self):
        x = np.random.standard_normal([1000, 10]) + 100.
        window = mlab.window_hanning
        window1 = mlab.window_hanning(np.ones(x.shape[0]))
        y, window2 = mlab.apply_window(x, window, axis=0, return_window=True)
        yt = mlab.apply_window(x, window1, axis=0, return_window=False)
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)
        assert_array_equal(window1, window2)

    def test_apply_window_hanning_2D_axis1(self):
        x = np.random.standard_normal([10, 1000]) + 100.
        window = mlab.window_hanning
        y = mlab.apply_window(x, window, axis=1, return_window=False)
        yt = np.zeros_like(x)
        for i in range(x.shape[0]):
            yt[i, :] = window(x[i, :])
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_hanning_2D__els1_axis1(self):
        x = np.random.standard_normal([10, 1000]) + 100.
        window = mlab.window_hanning(np.ones(x.shape[1]))
        window1 = mlab.window_hanning
        y = mlab.apply_window(x, window, axis=1, return_window=False)
        yt = np.zeros_like(x)
        for i in range(x.shape[0]):
            yt[i, :] = window1(x[i, :])
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_hanning_2D_els2_axis1(self):
        x = np.random.standard_normal([10, 1000]) + 100.
        window = mlab.window_hanning
        window1 = mlab.window_hanning(np.ones(x.shape[1]))
        y, window2 = mlab.apply_window(x, window, axis=1, return_window=True)
        yt = np.zeros_like(x)
        for i in range(x.shape[0]):
            yt[i, :] = window1 * x[i, :]
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)
        assert_array_equal(window1, window2)

    def test_apply_window_hanning_2D_els3_axis1(self):
        x = np.random.standard_normal([10, 1000]) + 100.
        window = mlab.window_hanning
        window1 = mlab.window_hanning(np.ones(x.shape[1]))
        y = mlab.apply_window(x, window, axis=1, return_window=False)
        yt = mlab.apply_window(x, window1, axis=1, return_window=False)
        assert yt.shape == y.shape
        assert x.shape == y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_stride_windows_hanning_2D_n13_noverlapn3_axis0(self):
        x = self.sig_rand
        window = mlab.window_hanning
        yi = mlab.stride_windows(x, n=13, noverlap=2, axis=0)
        y = mlab.apply_window(yi, window, axis=0, return_window=False)
        yt = self.check_window_apply_repeat(x, window, 13, 2)
        assert yt.shape == y.shape
        assert x.shape != y.shape
        assert_allclose(yt, y, atol=1e-06)

    def test_apply_window_hanning_2D_stack_axis1(self):
        ydata = np.arange(32)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ycontrol1 = mlab.apply_window(ydata1, mlab.window_hanning)
        ycontrol2 = mlab.window_hanning(ydata2)
        ydata = np.vstack([ydata1, ydata2])
        ycontrol = np.vstack([ycontrol1, ycontrol2])
        ydata = np.tile(ydata, (20, 1))
        ycontrol = np.tile(ycontrol, (20, 1))
        result = mlab.apply_window(ydata, mlab.window_hanning, axis=1,
                                   return_window=False)
        assert_allclose(ycontrol, result, atol=1e-08)

    def test_apply_window_hanning_2D_stack_windows_axis1(self):
        ydata = np.arange(32)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ycontrol1 = mlab.apply_window(ydata1, mlab.window_hanning)
        ycontrol2 = mlab.window_hanning(ydata2)
        ydata = np.vstack([ydata1, ydata2])
        ycontrol = np.vstack([ycontrol1, ycontrol2])
        ydata = np.tile(ydata, (20, 1))
        ycontrol = np.tile(ycontrol, (20, 1))
        result = mlab.apply_window(ydata, mlab.window_hanning, axis=1,
                                   return_window=False)
        assert_allclose(ycontrol, result, atol=1e-08)

    def test_apply_window_hanning_2D_stack_windows_axis1_unflatten(self):
        n = 32
        ydata = np.arange(n)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ycontrol1 = mlab.apply_window(ydata1, mlab.window_hanning)
        ycontrol2 = mlab.window_hanning(ydata2)
        ydata = np.vstack([ydata1, ydata2])
        ycontrol = np.vstack([ycontrol1, ycontrol2])
        ydata = np.tile(ydata, (20, 1))
        ycontrol = np.tile(ycontrol, (20, 1))
        ydata = ydata.flatten()
        ydata1 = mlab.stride_windows(ydata, 32, noverlap=0, axis=0)
        result = mlab.apply_window(ydata1, mlab.window_hanning, axis=0,
                                   return_window=False)
        assert_allclose(ycontrol.T, result, atol=1e-08)


class TestDetrend(object):
    def setup(self):
        np.random.seed(0)
        n = 1000
        x = np.linspace(0., 100, n)

        self.sig_zeros = np.zeros(n)

        self.sig_off = self.sig_zeros + 100.
        self.sig_slope = np.linspace(-10., 90., n)

        self.sig_slope_mean = x - x.mean()

        sig_rand = np.random.standard_normal(n)
        sig_sin = np.sin(x*2*np.pi/(n/100))

        sig_rand -= sig_rand.mean()
        sig_sin -= sig_sin.mean()

        self.sig_base = sig_rand + sig_sin

        self.atol = 1e-08

    def test_detrend_none_0D_zeros(self):
        input = 0.
        targ = input
        res = mlab.detrend_none(input)
        assert input == targ

    def test_detrend_none_0D_zeros_axis1(self):
        input = 0.
        targ = input
        res = mlab.detrend_none(input, axis=1)
        assert input == targ

    def test_detrend_str_none_0D_zeros(self):
        input = 0.
        targ = input
        res = mlab.detrend(input, key='none')
        assert input == targ

    def test_detrend_detrend_none_0D_zeros(self):
        input = 0.
        targ = input
        res = mlab.detrend(input, key=mlab.detrend_none)
        assert input == targ

    def test_detrend_none_0D_off(self):
        input = 5.5
        targ = input
        res = mlab.detrend_none(input)
        assert input == targ

    def test_detrend_none_1D_off(self):
        input = self.sig_off
        targ = input
        res = mlab.detrend_none(input)
        assert_array_equal(res, targ)

    def test_detrend_none_1D_slope(self):
        input = self.sig_slope
        targ = input
        res = mlab.detrend_none(input)
        assert_array_equal(res, targ)

    def test_detrend_none_1D_base(self):
        input = self.sig_base
        targ = input
        res = mlab.detrend_none(input)
        assert_array_equal(res, targ)

    def test_detrend_none_1D_base_slope_off_list(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = input.tolist()
        res = mlab.detrend_none(input.tolist())
        assert res == targ

    def test_detrend_none_2D(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        input = np.vstack(arri)
        targ = input
        res = mlab.detrend_none(input)
        assert_array_equal(res, targ)

    def test_detrend_none_2D_T(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        input = np.vstack(arri)
        targ = input
        res = mlab.detrend_none(input.T)
        assert_array_equal(res.T, targ)

    def test_detrend_mean_0D_zeros(self):
        input = 0.
        targ = 0.
        res = mlab.detrend_mean(input)
        assert_almost_equal(res, targ)

    def test_detrend_str_mean_0D_zeros(self):
        input = 0.
        targ = 0.
        res = mlab.detrend(input, key='mean')
        assert_almost_equal(res, targ)

    def test_detrend_detrend_mean_0D_zeros(self):
        input = 0.
        targ = 0.
        res = mlab.detrend(input, key=mlab.detrend_mean)
        assert_almost_equal(res, targ)

    def test_detrend_mean_0D_off(self):
        input = 5.5
        targ = 0.
        res = mlab.detrend_mean(input)
        assert_almost_equal(res, targ)

    def test_detrend_str_mean_0D_off(self):
        input = 5.5
        targ = 0.
        res = mlab.detrend(input, key='mean')
        assert_almost_equal(res, targ)

    def test_detrend_detrend_mean_0D_off(self):
        input = 5.5
        targ = 0.
        res = mlab.detrend(input, key=mlab.detrend_mean)
        assert_almost_equal(res, targ)

    def test_detrend_mean_1D_zeros(self):
        input = self.sig_zeros
        targ = self.sig_zeros
        res = mlab.detrend_mean(input)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_mean_1D_base(self):
        input = self.sig_base
        targ = self.sig_base
        res = mlab.detrend_mean(input)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_mean_1D_base_off(self):
        input = self.sig_base + self.sig_off
        targ = self.sig_base
        res = mlab.detrend_mean(input)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_mean_1D_base_slope(self):
        input = self.sig_base + self.sig_slope
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.detrend_mean(input)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_mean_1D_base_slope_off(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.detrend_mean(input)
        assert_allclose(res, targ, atol=1e-08)

    def test_detrend_mean_1D_base_slope_off_axis0(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.detrend_mean(input, axis=0)
        assert_allclose(res, targ, atol=1e-08)

    def test_detrend_mean_1D_base_slope_off_list(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.detrend_mean(input.tolist())
        assert_allclose(res, targ, atol=1e-08)

    def test_detrend_mean_1D_base_slope_off_list_axis0(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.detrend_mean(input.tolist(), axis=0)
        assert_allclose(res, targ, atol=1e-08)

    def test_demean_0D_off(self):
        input = 5.5
        targ = 0.
        res = mlab.demean(input, axis=None)
        assert_almost_equal(res, targ)

    def test_demean_1D_base_slope_off(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.demean(input)
        assert_allclose(res, targ, atol=1e-08)

    def test_demean_1D_base_slope_off_axis0(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.demean(input, axis=0)
        assert_allclose(res, targ, atol=1e-08)

    def test_demean_1D_base_slope_off_list(self):
        input = self.sig_base + self.sig_slope + self.sig_off
        targ = self.sig_base + self.sig_slope_mean
        res = mlab.demean(input.tolist())
        assert_allclose(res, targ, atol=1e-08)

    def test_detrend_mean_2D_default(self):
        arri = [self.sig_off,
                self.sig_base + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_base]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend_mean(input)
        assert_allclose(res, targ, atol=1e-08)

    def test_detrend_mean_2D_none(self):
        arri = [self.sig_off,
                self.sig_base + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_base]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend_mean(input, axis=None)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_detrend_mean_2D_none_T(self):
        arri = [self.sig_off,
                self.sig_base + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_base]
        input = np.vstack(arri).T
        targ = np.vstack(arrt)
        res = mlab.detrend_mean(input, axis=None)
        assert_allclose(res.T, targ,
                        atol=1e-08)

    def test_detrend_mean_2D_axis0(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri).T
        targ = np.vstack(arrt).T
        res = mlab.detrend_mean(input, axis=0)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_detrend_mean_2D_axis1(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend_mean(input, axis=1)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_detrend_mean_2D_axism1(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend_mean(input, axis=-1)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_detrend_2D_default(self):
        arri = [self.sig_off,
                self.sig_base + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_base]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend(input)
        assert_allclose(res, targ, atol=1e-08)

    def test_detrend_2D_none(self):
        arri = [self.sig_off,
                self.sig_base + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_base]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend(input, axis=None)
        assert_allclose(res, targ, atol=1e-08)

    def test_detrend_str_mean_2D_axis0(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri).T
        targ = np.vstack(arrt).T
        res = mlab.detrend(input, key='mean', axis=0)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_detrend_str_constant_2D_none_T(self):
        arri = [self.sig_off,
                self.sig_base + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_base]
        input = np.vstack(arri).T
        targ = np.vstack(arrt)
        res = mlab.detrend(input, key='constant', axis=None)
        assert_allclose(res.T, targ,
                        atol=1e-08)

    def test_detrend_str_default_2D_axis1(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend(input, key='default', axis=1)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_detrend_detrend_mean_2D_axis0(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri).T
        targ = np.vstack(arrt).T
        res = mlab.detrend(input, key=mlab.detrend_mean, axis=0)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_demean_2D_default(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri).T
        targ = np.vstack(arrt).T
        res = mlab.demean(input)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_demean_2D_none(self):
        arri = [self.sig_off,
                self.sig_base + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_base]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.demean(input, axis=None)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_demean_2D_axis0(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri).T
        targ = np.vstack(arrt).T
        res = mlab.demean(input, axis=0)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_demean_2D_axis1(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.demean(input, axis=1)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_demean_2D_axism1(self):
        arri = [self.sig_base,
                self.sig_base + self.sig_off,
                self.sig_base + self.sig_slope,
                self.sig_base + self.sig_off + self.sig_slope]
        arrt = [self.sig_base,
                self.sig_base,
                self.sig_base + self.sig_slope_mean,
                self.sig_base + self.sig_slope_mean]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.demean(input, axis=-1)
        assert_allclose(res, targ,
                        atol=1e-08)

    def test_detrend_bad_key_str_ValueError(self):
        input = self.sig_slope[np.newaxis]
        with pytest.raises(ValueError):
            mlab.detrend(input, key='spam')

    def test_detrend_bad_key_var_ValueError(self):
        input = self.sig_slope[np.newaxis]
        with pytest.raises(ValueError):
            mlab.detrend(input, key=5)

    def test_detrend_mean_0D_d0_ValueError(self):
        input = 5.5
        with pytest.raises(ValueError):
            mlab.detrend_mean(input, axis=0)

    def test_detrend_0D_d0_ValueError(self):
        input = 5.5
        with pytest.raises(ValueError):
            mlab.detrend(input, axis=0)

    def test_detrend_mean_1D_d1_ValueError(self):
        input = self.sig_slope
        with pytest.raises(ValueError):
            mlab.detrend_mean(input, axis=1)

    def test_detrend_1D_d1_ValueError(self):
        input = self.sig_slope
        with pytest.raises(ValueError):
            mlab.detrend(input, axis=1)

    def test_demean_1D_d1_ValueError(self):
        input = self.sig_slope
        with pytest.raises(ValueError):
            mlab.demean(input, axis=1)

    def test_detrend_mean_2D_d2_ValueError(self):
        input = self.sig_slope[np.newaxis]
        with pytest.raises(ValueError):
            mlab.detrend_mean(input, axis=2)

    def test_detrend_2D_d2_ValueError(self):
        input = self.sig_slope[np.newaxis]
        with pytest.raises(ValueError):
            mlab.detrend(input, axis=2)

    def test_demean_2D_d2_ValueError(self):
        input = self.sig_slope[np.newaxis]
        with pytest.raises(ValueError):
            mlab.demean(input, axis=2)

    def test_detrend_linear_0D_zeros(self):
        input = 0.
        targ = 0.
        res = mlab.detrend_linear(input)
        assert_almost_equal(res, targ)

    def test_detrend_linear_0D_off(self):
        input = 5.5
        targ = 0.
        res = mlab.detrend_linear(input)
        assert_almost_equal(res, targ)

    def test_detrend_str_linear_0D_off(self):
        input = 5.5
        targ = 0.
        res = mlab.detrend(input, key='linear')
        assert_almost_equal(res, targ)

    def test_detrend_detrend_linear_0D_off(self):
        input = 5.5
        targ = 0.
        res = mlab.detrend(input, key=mlab.detrend_linear)
        assert_almost_equal(res, targ)

    def test_detrend_linear_1d_off(self):
        input = self.sig_off
        targ = self.sig_zeros
        res = mlab.detrend_linear(input)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_linear_1d_slope(self):
        input = self.sig_slope
        targ = self.sig_zeros
        res = mlab.detrend_linear(input)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_linear_1d_slope_off(self):
        input = self.sig_slope + self.sig_off
        targ = self.sig_zeros
        res = mlab.detrend_linear(input)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_str_linear_1d_slope_off(self):
        input = self.sig_slope + self.sig_off
        targ = self.sig_zeros
        res = mlab.detrend(input, key='linear')
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_detrend_linear_1d_slope_off(self):
        input = self.sig_slope + self.sig_off
        targ = self.sig_zeros
        res = mlab.detrend(input, key=mlab.detrend_linear)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_linear_1d_slope_off_list(self):
        input = self.sig_slope + self.sig_off
        targ = self.sig_zeros
        res = mlab.detrend_linear(input.tolist())
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_linear_2D_ValueError(self):
        input = self.sig_slope[np.newaxis]
        with pytest.raises(ValueError):
            mlab.detrend_linear(input)

    def test_detrend_str_linear_2d_slope_off_axis0(self):
        arri = [self.sig_off,
                self.sig_slope,
                self.sig_slope + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_zeros,
                self.sig_zeros]
        input = np.vstack(arri).T
        targ = np.vstack(arrt).T
        res = mlab.detrend(input, key='linear', axis=0)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_detrend_linear_1d_slope_off_axis1(self):
        arri = [self.sig_off,
                self.sig_slope,
                self.sig_slope + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_zeros,
                self.sig_zeros]
        input = np.vstack(arri).T
        targ = np.vstack(arrt).T
        res = mlab.detrend(input, key=mlab.detrend_linear, axis=0)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_str_linear_2d_slope_off_axis0(self):
        arri = [self.sig_off,
                self.sig_slope,
                self.sig_slope + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_zeros,
                self.sig_zeros]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend(input, key='linear', axis=1)
        assert_allclose(res, targ, atol=self.atol)

    def test_detrend_detrend_linear_1d_slope_off_axis1(self):
        arri = [self.sig_off,
                self.sig_slope,
                self.sig_slope + self.sig_off]
        arrt = [self.sig_zeros,
                self.sig_zeros,
                self.sig_zeros]
        input = np.vstack(arri)
        targ = np.vstack(arrt)
        res = mlab.detrend(input, key=mlab.detrend_linear, axis=1)
        assert_allclose(res, targ, atol=self.atol)


@pytest.mark.parametrize('iscomplex', [False, True],
                         ids=['real', 'complex'], scope='class')
@pytest.mark.parametrize('sides', ['onesided', 'twosided', 'default'],
                         scope='class')
@pytest.mark.parametrize(
    'fstims,len_x,NFFT_density,nover_density,pad_to_density,pad_to_spectrum',
    [
        ([], None, -1, -1, -1, -1),
        ([4], None, -1, -1, -1, -1),
        ([4, 5, 10], None, -1, -1, -1, -1),
        ([], None, None, -1, -1, None),
        ([], None, -1, -1, None, None),
        ([], None, None, -1, None, None),
        ([], 1024, 512, -1, -1, 128),
        ([], 256, -1, -1, 33, 257),
        ([], 255, 33, -1, -1, None),
        ([], 256, 128, -1, 256, 256),
        ([], None, -1, 32, -1, -1),
   ],
    ids=[
        'nosig',
        'Fs4',
        'FsAll',
        'nosig_noNFFT',
        'nosig_nopad_to',
        'nosig_noNFFT_no_pad_to',
        'nosig_trim',
        'nosig_odd',
        'nosig_oddlen',
        'nosig_stretch',
        'nosig_overlap',
    ],
    scope='class')
class TestSpectral(object):
    @pytest.fixture(scope='class', autouse=True)
    def stim(self, request, fstims, iscomplex, sides, len_x, NFFT_density,
             nover_density, pad_to_density, pad_to_spectrum):
        Fs = 100.

        x = np.arange(0, 10, 1 / Fs)
        if len_x is not None:
            x = x[:len_x]

        # get the stimulus frequencies, defaulting to None
        fstims = [Fs / fstim for fstim in fstims]

        # get the constants, default to calculated values
        if NFFT_density is None:
            NFFT_density_real = 256
        elif NFFT_density < 0:
            NFFT_density_real = NFFT_density = 100
        else:
            NFFT_density_real = NFFT_density

        if nover_density is None:
            nover_density_real = 0
        elif nover_density < 0:
            nover_density_real = nover_density = NFFT_density_real // 2
        else:
            nover_density_real = nover_density

        if pad_to_density is None:
            pad_to_density_real = NFFT_density_real
        elif pad_to_density < 0:
            pad_to_density = int(2**np.ceil(np.log2(NFFT_density_real)))
            pad_to_density_real = pad_to_density
        else:
            pad_to_density_real = pad_to_density

        if pad_to_spectrum is None:
            pad_to_spectrum_real = len(x)
        elif pad_to_spectrum < 0:
            pad_to_spectrum_real = pad_to_spectrum = len(x)
        else:
            pad_to_spectrum_real = pad_to_spectrum

        if pad_to_spectrum is None:
            NFFT_spectrum_real = NFFT_spectrum = pad_to_spectrum_real
        else:
            NFFT_spectrum_real = NFFT_spectrum = len(x)
        nover_spectrum_real = nover_spectrum = 0

        NFFT_specgram = NFFT_density
        nover_specgram = nover_density
        pad_to_specgram = pad_to_density
        NFFT_specgram_real = NFFT_density_real
        nover_specgram_real = nover_density_real

        if sides == 'onesided' or (sides == 'default' and not iscomplex):
            # frequencies for specgram, psd, and csd
            # need to handle even and odd differently
            if pad_to_density_real % 2:
                freqs_density = np.linspace(0, Fs / 2,
                                            num=pad_to_density_real,
                                            endpoint=False)[::2]
            else:
                freqs_density = np.linspace(0, Fs / 2,
                                            num=pad_to_density_real // 2 + 1)

            # frequencies for complex, magnitude, angle, and phase spectrums
            # need to handle even and odd differently
            if pad_to_spectrum_real % 2:
                freqs_spectrum = np.linspace(0, Fs / 2,
                                             num=pad_to_spectrum_real,
                                             endpoint=False)[::2]
            else:
                freqs_spectrum = np.linspace(0, Fs / 2,
                                             num=pad_to_spectrum_real // 2 + 1)
        else:
            # frequencies for specgram, psd, and csd
            # need to handle even and odd differentl
            if pad_to_density_real % 2:
                freqs_density = np.linspace(-Fs / 2, Fs / 2,
                                            num=2 * pad_to_density_real,
                                            endpoint=False)[1::2]
            else:
                freqs_density = np.linspace(-Fs / 2, Fs / 2,
                                            num=pad_to_density_real,
                                            endpoint=False)

            # frequencies for complex, magnitude, angle, and phase spectrums
            # need to handle even and odd differently
            if pad_to_spectrum_real % 2:
                freqs_spectrum = np.linspace(-Fs / 2, Fs / 2,
                                             num=2 * pad_to_spectrum_real,
                                             endpoint=False)[1::2]
            else:
                freqs_spectrum = np.linspace(-Fs / 2, Fs / 2,
                                             num=pad_to_spectrum_real,
                                             endpoint=False)

        freqs_specgram = freqs_density
        # time points for specgram
        t_start = NFFT_specgram_real // 2
        t_stop = len(x) - NFFT_specgram_real // 2 + 1
        t_step = NFFT_specgram_real - nover_specgram_real
        t_specgram = x[t_start:t_stop:t_step]
        if NFFT_specgram_real % 2:
            t_specgram += 1 / Fs / 2
        if len(t_specgram) == 0:
            t_specgram = np.array([NFFT_specgram_real / (2 * Fs)])
        t_spectrum = np.array([NFFT_spectrum_real / (2 * Fs)])
        t_density = t_specgram

        y = np.zeros_like(x)
        for i, fstim in enumerate(fstims):
            y += np.sin(fstim * x * np.pi * 2) * 10**i

        if iscomplex:
            y = y.astype('complex')

        # Interestingly, the instance on which this fixture is called is not
        # the same as the one on which a test is run. So we need to modify the
        # class itself when using a class-scoped fixture.
        cls = request.cls

        cls.Fs = Fs
        cls.sides = sides
        cls.fstims = fstims

        cls.NFFT_density = NFFT_density
        cls.nover_density = nover_density
        cls.pad_to_density = pad_to_density

        cls.NFFT_spectrum = NFFT_spectrum
        cls.nover_spectrum = nover_spectrum
        cls.pad_to_spectrum = pad_to_spectrum

        cls.NFFT_specgram = NFFT_specgram
        cls.nover_specgram = nover_specgram
        cls.pad_to_specgram = pad_to_specgram

        cls.t_specgram = t_specgram
        cls.t_density = t_density
        cls.t_spectrum = t_spectrum
        cls.y = y

        cls.freqs_density = freqs_density
        cls.freqs_spectrum = freqs_spectrum
        cls.freqs_specgram = freqs_specgram

        cls.NFFT_density_real = NFFT_density_real

    def check_freqs(self, vals, targfreqs, resfreqs, fstims):
        assert resfreqs.argmin() == 0
        assert resfreqs.argmax() == len(resfreqs)-1
        assert_allclose(resfreqs, targfreqs, atol=1e-06)
        for fstim in fstims:
            i = np.abs(resfreqs - fstim).argmin()
            assert vals[i] > vals[i+2]
            assert vals[i] > vals[i-2]

    def check_maxfreq(self, spec, fsp, fstims):
        # skip the test if there are no frequencies
        if len(fstims) == 0:
            return

        # if twosided, do the test for each side
        if fsp.min() < 0:
            fspa = np.abs(fsp)
            zeroind = fspa.argmin()
            self.check_maxfreq(spec[:zeroind], fspa[:zeroind], fstims)
            self.check_maxfreq(spec[zeroind:], fspa[zeroind:], fstims)
            return

        fstimst = fstims[:]
        spect = spec.copy()

        # go through each peak and make sure it is correctly the maximum peak
        while fstimst:
            maxind = spect.argmax()
            maxfreq = fsp[maxind]
            assert_almost_equal(maxfreq, fstimst[-1])
            del fstimst[-1]
            spect[maxind-5:maxind+5] = 0

    def test_spectral_helper_raises_complex_same_data(self):
        # test that mode 'complex' cannot be used if x is not y
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, y=self.y+1, mode='complex')

    def test_spectral_helper_raises_magnitude_same_data(self):
        # test that mode 'magnitude' cannot be used if x is not y
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, y=self.y+1, mode='magnitude')

    def test_spectral_helper_raises_angle_same_data(self):
        # test that mode 'angle' cannot be used if x is not y
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, y=self.y+1, mode='angle')

    def test_spectral_helper_raises_phase_same_data(self):
        # test that mode 'phase' cannot be used if x is not y
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, y=self.y+1, mode='phase')

    def test_spectral_helper_raises_unknown_mode(self):
        # test that unknown value for mode cannot be used
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, mode='spam')

    def test_spectral_helper_raises_unknown_sides(self):
        # test that unknown value for sides cannot be used
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, y=self.y, sides='eggs')

    def test_spectral_helper_raises_noverlap_gt_NFFT(self):
        # test that noverlap cannot be larger than NFFT
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, y=self.y, NFFT=10, noverlap=20)

    def test_spectral_helper_raises_noverlap_eq_NFFT(self):
        # test that noverlap cannot be equal to NFFT
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, NFFT=10, noverlap=10)

    def test_spectral_helper_raises_winlen_ne_NFFT(self):
        # test that the window length cannot be different from NFFT
        with pytest.raises(ValueError):
            mlab._spectral_helper(x=self.y, y=self.y, NFFT=10,
                                  window=np.ones(9))

    def test_single_spectrum_helper_raises_mode_default(self):
        # test that mode 'default' cannot be used with _single_spectrum_helper
        with pytest.raises(ValueError):
            mlab._single_spectrum_helper(x=self.y, mode='default')

    def test_single_spectrum_helper_raises_mode_psd(self):
        # test that mode 'psd' cannot be used with _single_spectrum_helper
        with pytest.raises(ValueError):
            mlab._single_spectrum_helper(x=self.y, mode='psd')

    def test_spectral_helper_psd(self):
        freqs = self.freqs_density
        spec, fsp, t = mlab._spectral_helper(x=self.y, y=self.y,
                                             NFFT=self.NFFT_density,
                                             Fs=self.Fs,
                                             noverlap=self.nover_density,
                                             pad_to=self.pad_to_density,
                                             sides=self.sides,
                                             mode='psd')

        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_density, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]

    def test_spectral_helper_magnitude_specgram(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab._spectral_helper(x=self.y, y=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='magnitude')

        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]

    def test_spectral_helper_magnitude_magnitude_spectrum(self):
        freqs = self.freqs_spectrum
        spec, fsp, t = mlab._spectral_helper(x=self.y, y=self.y,
                                             NFFT=self.NFFT_spectrum,
                                             Fs=self.Fs,
                                             noverlap=self.nover_spectrum,
                                             pad_to=self.pad_to_spectrum,
                                             sides=self.sides,
                                             mode='magnitude')

        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_spectrum, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == 1

    def test_csd(self):
        freqs = self.freqs_density
        spec, fsp = mlab.csd(x=self.y, y=self.y+1,
                             NFFT=self.NFFT_density,
                             Fs=self.Fs,
                             noverlap=self.nover_density,
                             pad_to=self.pad_to_density,
                             sides=self.sides)
        assert_allclose(fsp, freqs, atol=1e-06)
        assert spec.shape == freqs.shape

    def test_csd_padding(self):
        """Test zero padding of csd(). """
        if self.NFFT_density is None:  # for derived classes
            return
        sargs = dict(x=self.y, y=self.y+1, Fs=self.Fs, window=mlab.window_none,
                     sides=self.sides)

        spec0, _ = mlab.csd(NFFT=self.NFFT_density, **sargs)
        spec1, _ = mlab.csd(NFFT=self.NFFT_density*2, **sargs)
        assert_almost_equal(np.sum(np.conjugate(spec0)*spec0).real,
                            np.sum(np.conjugate(spec1/2)*spec1/2).real)

    def test_psd(self):
        freqs = self.freqs_density
        spec, fsp = mlab.psd(x=self.y,
                             NFFT=self.NFFT_density,
                             Fs=self.Fs,
                             noverlap=self.nover_density,
                             pad_to=self.pad_to_density,
                             sides=self.sides)
        assert spec.shape == freqs.shape
        self.check_freqs(spec, freqs, fsp, self.fstims)

    def test_psd_detrend_mean_func_offset(self):
        if self.NFFT_density is None:
            return
        freqs = self.freqs_density
        ydata = np.zeros(self.NFFT_density)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ydata = np.vstack([ydata1, ydata2])
        ydata = np.tile(ydata, (20, 1))
        ydatab = ydata.T.flatten()
        ydata = ydata.flatten()
        ycontrol = np.zeros_like(ydata)
        spec_g, fsp_g = mlab.psd(x=ydata,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend=mlab.detrend_mean)
        spec_b, fsp_b = mlab.psd(x=ydatab,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend=mlab.detrend_mean)
        spec_c, fsp_c = mlab.psd(x=ycontrol,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides)
        assert_array_equal(fsp_g, fsp_c)
        assert_array_equal(fsp_b, fsp_c)
        assert_allclose(spec_g, spec_c, atol=1e-08)
        # these should not be almost equal
        with pytest.raises(AssertionError):
            assert_allclose(spec_b, spec_c, atol=1e-08)

    def test_psd_detrend_mean_str_offset(self):
        if self.NFFT_density is None:
            return
        freqs = self.freqs_density
        ydata = np.zeros(self.NFFT_density)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ydata = np.vstack([ydata1, ydata2])
        ydata = np.tile(ydata, (20, 1))
        ydatab = ydata.T.flatten()
        ydata = ydata.flatten()
        ycontrol = np.zeros_like(ydata)
        spec_g, fsp_g = mlab.psd(x=ydata,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend='mean')
        spec_b, fsp_b = mlab.psd(x=ydatab,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend='mean')
        spec_c, fsp_c = mlab.psd(x=ycontrol,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides)
        assert_array_equal(fsp_g, fsp_c)
        assert_array_equal(fsp_b, fsp_c)
        assert_allclose(spec_g, spec_c, atol=1e-08)
        # these should not be almost equal
        with pytest.raises(AssertionError):
            assert_allclose(spec_b, spec_c, atol=1e-08)

    def test_psd_detrend_linear_func_trend(self):
        if self.NFFT_density is None:
            return
        freqs = self.freqs_density
        ydata = np.arange(self.NFFT_density)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ydata = np.vstack([ydata1, ydata2])
        ydata = np.tile(ydata, (20, 1))
        ydatab = ydata.T.flatten()
        ydata = ydata.flatten()
        ycontrol = np.zeros_like(ydata)
        spec_g, fsp_g = mlab.psd(x=ydata,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend=mlab.detrend_linear)
        spec_b, fsp_b = mlab.psd(x=ydatab,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend=mlab.detrend_linear)
        spec_c, fsp_c = mlab.psd(x=ycontrol,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides)
        assert_array_equal(fsp_g, fsp_c)
        assert_array_equal(fsp_b, fsp_c)
        assert_allclose(spec_g, spec_c, atol=1e-08)
        # these should not be almost equal
        with pytest.raises(AssertionError):
            assert_allclose(spec_b, spec_c, atol=1e-08)

    def test_psd_detrend_linear_str_trend(self):
        if self.NFFT_density is None:
            return
        freqs = self.freqs_density
        ydata = np.arange(self.NFFT_density)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ydata = np.vstack([ydata1, ydata2])
        ydata = np.tile(ydata, (20, 1))
        ydatab = ydata.T.flatten()
        ydata = ydata.flatten()
        ycontrol = np.zeros_like(ydata)
        spec_g, fsp_g = mlab.psd(x=ydata,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend='linear')
        spec_b, fsp_b = mlab.psd(x=ydatab,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend='linear')
        spec_c, fsp_c = mlab.psd(x=ycontrol,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides)
        assert_array_equal(fsp_g, fsp_c)
        assert_array_equal(fsp_b, fsp_c)
        assert_allclose(spec_g, spec_c, atol=1e-08)
        # these should not be almost equal
        with pytest.raises(AssertionError):
            assert_allclose(spec_b, spec_c, atol=1e-08)

    def test_psd_window_hanning(self):
        if self.NFFT_density is None:
            return
        freqs = self.freqs_density
        ydata = np.arange(self.NFFT_density)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ycontrol1, windowVals = mlab.apply_window(ydata1,
                                                  mlab.window_hanning,
                                                  return_window=True)
        ycontrol2 = mlab.window_hanning(ydata2)
        ydata = np.vstack([ydata1, ydata2])
        ycontrol = np.vstack([ycontrol1, ycontrol2])
        ydata = np.tile(ydata, (20, 1))
        ycontrol = np.tile(ycontrol, (20, 1))
        ydatab = ydata.T.flatten()
        ydataf = ydata.flatten()
        ycontrol = ycontrol.flatten()
        spec_g, fsp_g = mlab.psd(x=ydataf,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 window=mlab.window_hanning)
        spec_b, fsp_b = mlab.psd(x=ydatab,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 window=mlab.window_hanning)
        spec_c, fsp_c = mlab.psd(x=ycontrol,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 window=mlab.window_none)
        spec_c *= len(ycontrol1)/(np.abs(windowVals)**2).sum()
        assert_array_equal(fsp_g, fsp_c)
        assert_array_equal(fsp_b, fsp_c)
        assert_allclose(spec_g, spec_c, atol=1e-08)
        # these should not be almost equal
        with pytest.raises(AssertionError):
            assert_allclose(spec_b, spec_c, atol=1e-08)

    def test_psd_window_hanning_detrend_linear(self):
        if self.NFFT_density is None:
            return
        freqs = self.freqs_density
        ydata = np.arange(self.NFFT_density)
        ycontrol = np.zeros(self.NFFT_density)
        ydata1 = ydata+5
        ydata2 = ydata+3.3
        ycontrol1 = ycontrol
        ycontrol2 = ycontrol
        ycontrol1, windowVals = mlab.apply_window(ycontrol1,
                                                  mlab.window_hanning,
                                                  return_window=True)
        ycontrol2 = mlab.window_hanning(ycontrol2)
        ydata = np.vstack([ydata1, ydata2])
        ycontrol = np.vstack([ycontrol1, ycontrol2])
        ydata = np.tile(ydata, (20, 1))
        ycontrol = np.tile(ycontrol, (20, 1))
        ydatab = ydata.T.flatten()
        ydataf = ydata.flatten()
        ycontrol = ycontrol.flatten()
        spec_g, fsp_g = mlab.psd(x=ydataf,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend=mlab.detrend_linear,
                                 window=mlab.window_hanning)
        spec_b, fsp_b = mlab.psd(x=ydatab,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 detrend=mlab.detrend_linear,
                                 window=mlab.window_hanning)
        spec_c, fsp_c = mlab.psd(x=ycontrol,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=0,
                                 sides=self.sides,
                                 window=mlab.window_none)
        spec_c *= len(ycontrol1)/(np.abs(windowVals)**2).sum()
        assert_array_equal(fsp_g, fsp_c)
        assert_array_equal(fsp_b, fsp_c)
        assert_allclose(spec_g, spec_c, atol=1e-08)
        # these should not be almost equal
        with pytest.raises(AssertionError):
            assert_allclose(spec_b, spec_c, atol=1e-08)

    def test_psd_windowarray(self):
        freqs = self.freqs_density
        spec, fsp = mlab.psd(x=self.y,
                             NFFT=self.NFFT_density,
                             Fs=self.Fs,
                             noverlap=self.nover_density,
                             pad_to=self.pad_to_density,
                             sides=self.sides,
                             window=np.ones(self.NFFT_density_real))
        assert_allclose(fsp, freqs, atol=1e-06)
        assert spec.shape == freqs.shape

    def test_psd_windowarray_scale_by_freq(self):
        freqs = self.freqs_density
        win = mlab.window_hanning(np.ones(self.NFFT_density_real))

        spec, fsp = mlab.psd(x=self.y,
                             NFFT=self.NFFT_density,
                             Fs=self.Fs,
                             noverlap=self.nover_density,
                             pad_to=self.pad_to_density,
                             sides=self.sides,
                             window=mlab.window_hanning)
        spec_s, fsp_s = mlab.psd(x=self.y,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=self.nover_density,
                                 pad_to=self.pad_to_density,
                                 sides=self.sides,
                                 window=mlab.window_hanning,
                                 scale_by_freq=True)
        spec_n, fsp_n = mlab.psd(x=self.y,
                                 NFFT=self.NFFT_density,
                                 Fs=self.Fs,
                                 noverlap=self.nover_density,
                                 pad_to=self.pad_to_density,
                                 sides=self.sides,
                                 window=mlab.window_hanning,
                                 scale_by_freq=False)
        assert_array_equal(fsp, fsp_s)
        assert_array_equal(fsp, fsp_n)
        assert_array_equal(spec, spec_s)
        assert_allclose(spec_s*(win**2).sum(),
                        spec_n/self.Fs*win.sum()**2,
                        atol=1e-08)

    def test_complex_spectrum(self):
        freqs = self.freqs_spectrum
        spec, fsp = mlab.complex_spectrum(x=self.y,
                                          Fs=self.Fs,
                                          sides=self.sides,
                                          pad_to=self.pad_to_spectrum)
        assert_allclose(fsp, freqs, atol=1e-06)
        assert spec.shape == freqs.shape

    def test_magnitude_spectrum(self):
        freqs = self.freqs_spectrum
        spec, fsp = mlab.magnitude_spectrum(x=self.y,
                                            Fs=self.Fs,
                                            sides=self.sides,
                                            pad_to=self.pad_to_spectrum)
        assert spec.shape == freqs.shape
        self.check_maxfreq(spec, fsp, self.fstims)
        self.check_freqs(spec, freqs, fsp, self.fstims)

    def test_angle_spectrum(self):
        freqs = self.freqs_spectrum
        spec, fsp = mlab.angle_spectrum(x=self.y,
                                        Fs=self.Fs,
                                        sides=self.sides,
                                        pad_to=self.pad_to_spectrum)
        assert_allclose(fsp, freqs, atol=1e-06)
        assert spec.shape == freqs.shape

    def test_phase_spectrum(self):
        freqs = self.freqs_spectrum
        spec, fsp = mlab.phase_spectrum(x=self.y,
                                        Fs=self.Fs,
                                        sides=self.sides,
                                        pad_to=self.pad_to_spectrum)
        assert_allclose(fsp, freqs, atol=1e-06)
        assert spec.shape == freqs.shape

    def test_specgram_auto(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab.specgram(x=self.y,
                                     NFFT=self.NFFT_specgram,
                                     Fs=self.Fs,
                                     noverlap=self.nover_specgram,
                                     pad_to=self.pad_to_specgram,
                                     sides=self.sides)
        specm = np.mean(spec, axis=1)

        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]

        # since we are using a single freq, all time slices
        # should be about the same
        if np.abs(spec.max()) != 0:
            assert_allclose(np.diff(spec, axis=1).max()/np.abs(spec.max()), 0,
                            atol=1e-02)
        self.check_freqs(specm, freqs, fsp, self.fstims)

    def test_specgram_default(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab.specgram(x=self.y,
                                     NFFT=self.NFFT_specgram,
                                     Fs=self.Fs,
                                     noverlap=self.nover_specgram,
                                     pad_to=self.pad_to_specgram,
                                     sides=self.sides,
                                     mode='default')
        specm = np.mean(spec, axis=1)

        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]

        # since we are using a single freq, all time slices
        # should be about the same
        if np.abs(spec.max()) != 0:
            assert_allclose(np.diff(spec, axis=1).max()/np.abs(spec.max()), 0,
                            atol=1e-02)
        self.check_freqs(specm, freqs, fsp, self.fstims)

    def test_specgram_psd(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab.specgram(x=self.y,
                                     NFFT=self.NFFT_specgram,
                                     Fs=self.Fs,
                                     noverlap=self.nover_specgram,
                                     pad_to=self.pad_to_specgram,
                                     sides=self.sides,
                                     mode='psd')
        specm = np.mean(spec, axis=1)

        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]
        # since we are using a single freq, all time slices
        # should be about the same
        if np.abs(spec.max()) != 0:
            assert_allclose(np.diff(spec, axis=1).max()/np.abs(spec.max()), 0,
                            atol=1e-02)
        self.check_freqs(specm, freqs, fsp, self.fstims)

    def test_specgram_complex(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab.specgram(x=self.y,
                                     NFFT=self.NFFT_specgram,
                                     Fs=self.Fs,
                                     noverlap=self.nover_specgram,
                                     pad_to=self.pad_to_specgram,
                                     sides=self.sides,
                                     mode='complex')
        specm = np.mean(np.abs(spec), axis=1)
        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]

        self.check_freqs(specm, freqs, fsp, self.fstims)

    def test_specgram_magnitude(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab.specgram(x=self.y,
                                     NFFT=self.NFFT_specgram,
                                     Fs=self.Fs,
                                     noverlap=self.nover_specgram,
                                     pad_to=self.pad_to_specgram,
                                     sides=self.sides,
                                     mode='magnitude')
        specm = np.mean(spec, axis=1)
        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]
        # since we are using a single freq, all time slices
        # should be about the same
        if np.abs(spec.max()) != 0:
            assert_allclose(np.diff(spec, axis=1).max()/np.abs(spec.max()), 0,
                            atol=1e-02)
        self.check_freqs(specm, freqs, fsp, self.fstims)

    def test_specgram_angle(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab.specgram(x=self.y,
                                     NFFT=self.NFFT_specgram,
                                     Fs=self.Fs,
                                     noverlap=self.nover_specgram,
                                     pad_to=self.pad_to_specgram,
                                     sides=self.sides,
                                     mode='angle')
        specm = np.mean(spec, axis=1)
        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]

    def test_specgram_phase(self):
        freqs = self.freqs_specgram
        spec, fsp, t = mlab.specgram(x=self.y,
                                     NFFT=self.NFFT_specgram,
                                     Fs=self.Fs,
                                     noverlap=self.nover_specgram,
                                     pad_to=self.pad_to_specgram,
                                     sides=self.sides,
                                     mode='phase')
        specm = np.mean(spec, axis=1)

        assert_allclose(fsp, freqs, atol=1e-06)
        assert_allclose(t, self.t_specgram, atol=1e-06)

        assert spec.shape[0] == freqs.shape[0]
        assert spec.shape[1] == self.t_specgram.shape[0]

    def test_specgram_warn_only1seg(self):
        """Warning should be raised if len(x) <= NFFT. """
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", category=UserWarning)
            mlab.specgram(x=self.y, NFFT=len(self.y), Fs=self.Fs)
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert str(w[0].message).startswith("Only one segment is calculated")

    def test_psd_csd_equal(self):
        freqs = self.freqs_density
        Pxx, freqsxx = mlab.psd(x=self.y,
                                NFFT=self.NFFT_density,
                                Fs=self.Fs,
                                noverlap=self.nover_density,
                                pad_to=self.pad_to_density,
                                sides=self.sides)
        Pxy, freqsxy = mlab.csd(x=self.y, y=self.y,
                                NFFT=self.NFFT_density,
                                Fs=self.Fs,
                                noverlap=self.nover_density,
                                pad_to=self.pad_to_density,
                                sides=self.sides)
        assert_array_almost_equal_nulp(Pxx, Pxy)
        assert_array_equal(freqsxx, freqsxy)

    def test_specgram_auto_default_equal(self):
        '''test that mlab.specgram without mode and with mode 'default' and
        'psd' are all the same'''
        freqs = self.freqs_specgram
        speca, freqspeca, ta = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides)
        specb, freqspecb, tb = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='default')
        assert_array_equal(speca, specb)
        assert_array_equal(freqspeca, freqspecb)
        assert_array_equal(ta, tb)

    def test_specgram_auto_psd_equal(self):
        '''test that mlab.specgram without mode and with mode 'default' and
        'psd' are all the same'''
        freqs = self.freqs_specgram
        speca, freqspeca, ta = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides)
        specc, freqspecc, tc = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='psd')
        assert_array_equal(speca, specc)
        assert_array_equal(freqspeca, freqspecc)
        assert_array_equal(ta, tc)

    def test_specgram_complex_mag_equivalent(self):
        freqs = self.freqs_specgram
        specc, freqspecc, tc = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='complex')
        specm, freqspecm, tm = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='magnitude')

        assert_array_equal(freqspecc, freqspecm)
        assert_array_equal(tc, tm)
        assert_allclose(np.abs(specc), specm, atol=1e-06)

    def test_specgram_complex_angle_equivalent(self):
        freqs = self.freqs_specgram
        specc, freqspecc, tc = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='complex')
        speca, freqspeca, ta = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='angle')

        assert_array_equal(freqspecc, freqspeca)
        assert_array_equal(tc, ta)
        assert_allclose(np.angle(specc), speca, atol=1e-06)

    def test_specgram_complex_phase_equivalent(self):
        freqs = self.freqs_specgram
        specc, freqspecc, tc = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='complex')
        specp, freqspecp, tp = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='phase')

        assert_array_equal(freqspecc, freqspecp)
        assert_array_equal(tc, tp)
        assert_allclose(np.unwrap(np.angle(specc), axis=0), specp,
                        atol=1e-06)

    def test_specgram_angle_phase_equivalent(self):
        freqs = self.freqs_specgram
        speca, freqspeca, ta = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='angle')
        specp, freqspecp, tp = mlab.specgram(x=self.y,
                                             NFFT=self.NFFT_specgram,
                                             Fs=self.Fs,
                                             noverlap=self.nover_specgram,
                                             pad_to=self.pad_to_specgram,
                                             sides=self.sides,
                                             mode='phase')

        assert_array_equal(freqspeca, freqspecp)
        assert_array_equal(ta, tp)
        assert_allclose(np.unwrap(speca, axis=0), specp,
                        atol=1e-06)

    def test_psd_windowarray_equal(self):
        freqs = self.freqs_density
        win = mlab.window_hanning(np.ones(self.NFFT_density_real))
        speca, fspa = mlab.psd(x=self.y,
                               NFFT=self.NFFT_density,
                               Fs=self.Fs,
                               noverlap=self.nover_density,
                               pad_to=self.pad_to_density,
                               sides=self.sides,
                               window=win)
        specb, fspb = mlab.psd(x=self.y,
                               NFFT=self.NFFT_density,
                               Fs=self.Fs,
                               noverlap=self.nover_density,
                               pad_to=self.pad_to_density,
                               sides=self.sides)
        assert_array_equal(fspa, fspb)
        assert_allclose(speca, specb, atol=1e-08)


# extra test for cohere...
def test_cohere():
    N = 1024
    np.random.seed(19680801)
    x = np.random.randn(N)
    # phase offset
    y = np.roll(x, 20)
    # high-freq roll-off
    y = np.convolve(y, np.ones(20) / 20., mode='same')
    cohsq, f = mlab.cohere(x, y, NFFT=256, Fs=2, noverlap=128)
    assert_allclose(np.mean(cohsq), 0.837, atol=1.e-3)
    assert np.isreal(np.mean(cohsq))


def test_griddata_linear():
    # z is a linear function of x and y.
    def get_z(x, y):
        return 3.0*x - y

    # Passing 1D xi and yi arrays to griddata.
    x = np.asarray([0.0, 1.0, 0.0, 1.0, 0.5])
    y = np.asarray([0.0, 0.0, 1.0, 1.0, 0.5])
    z = get_z(x, y)
    xi = [0.2, 0.4, 0.6, 0.8]
    yi = [0.1, 0.3, 0.7, 0.9]
    with pytest.warns(MatplotlibDeprecationWarning):
        zi = mlab.griddata(x, y, z, xi, yi, interp='linear')
    xi, yi = np.meshgrid(xi, yi)
    np.testing.assert_array_almost_equal(zi, get_z(xi, yi))

    # Passing 2D xi and yi arrays to griddata.
    with pytest.warns(MatplotlibDeprecationWarning):
        zi = mlab.griddata(x, y, z, xi, yi, interp='linear')
    np.testing.assert_array_almost_equal(zi, get_z(xi, yi))

    # Masking z array.
    z_masked = np.ma.array(z, mask=[False, False, False, True, False])
    correct_zi_masked = np.ma.masked_where(xi + yi > 1.0, get_z(xi, yi))
    with pytest.warns(MatplotlibDeprecationWarning):
        zi = mlab.griddata(x, y, z_masked, xi, yi, interp='linear')
    matest.assert_array_almost_equal(zi, correct_zi_masked)
    np.testing.assert_array_equal(np.ma.getmask(zi),
                                  np.ma.getmask(correct_zi_masked))


@pytest.mark.xfail(not HAS_NATGRID, reason='natgrid not installed')
def test_griddata_nn():
    # z is a linear function of x and y.
    def get_z(x, y):
        return 3.0*x - y

    # Passing 1D xi and yi arrays to griddata.
    x = np.asarray([0.0, 1.0, 0.0, 1.0, 0.5])
    y = np.asarray([0.0, 0.0, 1.0, 1.0, 0.5])
    z = get_z(x, y)
    xi = [0.2, 0.4, 0.6, 0.8]
    yi = [0.1, 0.3, 0.7, 0.9]
    correct_zi = [[0.49999252, 1.0999978, 1.7000030, 2.3000080],
                  [0.29999208, 0.8999978, 1.5000029, 2.1000059],
                  [-0.1000099, 0.4999943, 1.0999964, 1.6999979],
                  [-0.3000128, 0.2999894, 0.8999913, 1.4999933]]
    with pytest.warns(MatplotlibDeprecationWarning):
        zi = mlab.griddata(x, y, z, xi, yi, interp='nn')
    np.testing.assert_array_almost_equal(zi, correct_zi, 5)

    with pytest.warns(MatplotlibDeprecationWarning):
        # Decreasing xi or yi should raise ValueError.
        with pytest.raises(ValueError):
            mlab.griddata(x, y, z, xi[::-1], yi, interp='nn')
        with pytest.raises(ValueError):
            mlab.griddata(x, y, z, xi, yi[::-1], interp='nn')

    # Passing 2D xi and yi arrays to griddata.
    xi, yi = np.meshgrid(xi, yi)
    with pytest.warns(MatplotlibDeprecationWarning):
        zi = mlab.griddata(x, y, z, xi, yi, interp='nn')
    np.testing.assert_array_almost_equal(zi, correct_zi, 5)

    # Masking z array.
    z_masked = np.ma.array(z, mask=[False, False, False, True, False])
    correct_zi_masked = np.ma.masked_where(xi + yi > 1.0, correct_zi)
    with pytest.warns(MatplotlibDeprecationWarning):
        zi = mlab.griddata(x, y, z_masked, xi, yi, interp='nn')
    np.testing.assert_array_almost_equal(zi, correct_zi_masked, 5)
    np.testing.assert_array_equal(np.ma.getmask(zi),
                                  np.ma.getmask(correct_zi_masked))


#*****************************************************************
# These Tests where taken from SCIPY with some minor modifications
# this can be retrieved from:
# https://github.com/scipy/scipy/blob/master/scipy/stats/tests/test_kdeoth.py
#*****************************************************************

class TestGaussianKDE(object):

    def test_kde_integer_input(self):
        """Regression test for #1181."""
        x1 = np.arange(5)
        kde = mlab.GaussianKDE(x1)
        y_expected = [0.13480721, 0.18222869, 0.19514935, 0.18222869,
                      0.13480721]
        np.testing.assert_array_almost_equal(kde(x1), y_expected, decimal=6)

    def test_gaussian_kde_covariance_caching(self):
        x1 = np.array([-7, -5, 1, 4, 5], dtype=float)
        xs = np.linspace(-10, 10, num=5)
        # These expected values are from scipy 0.10, before some changes to
        # gaussian_kde. They were not compared with any external reference.
        y_expected = [0.02463386, 0.04689208, 0.05395444, 0.05337754,
                      0.01664475]

        # set it to the default bandwidth.
        kde2 = mlab.GaussianKDE(x1, 'scott')
        y2 = kde2(xs)

        np.testing.assert_array_almost_equal(y_expected, y2, decimal=7)

    def test_kde_bandwidth_method(self):

        np.random.seed(8765678)
        n_basesample = 50
        xn = np.random.randn(n_basesample)

        # Default
        gkde = mlab.GaussianKDE(xn)
        # Supply a callable
        gkde2 = mlab.GaussianKDE(xn, 'scott')
        # Supply a scalar
        gkde3 = mlab.GaussianKDE(xn, bw_method=gkde.factor)

        xs = np.linspace(-7, 7, 51)
        kdepdf = gkde.evaluate(xs)
        kdepdf2 = gkde2.evaluate(xs)
        assert kdepdf.all() == kdepdf2.all()
        kdepdf3 = gkde3.evaluate(xs)
        assert kdepdf.all() == kdepdf3.all()


class TestGaussianKDECustom(object):
    def test_no_data(self):
        """Pass no data into the GaussianKDE class."""
        with pytest.raises(ValueError):
            mlab.GaussianKDE([])

    def test_single_dataset_element(self):
        """Pass a single dataset element into the GaussianKDE class."""
        with pytest.raises(ValueError):
            mlab.GaussianKDE([42])

    def test_silverman_multidim_dataset(self):
        """Use a multi-dimensional array as the dataset and test silverman's
        output"""
        x1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        with pytest.raises(np.linalg.LinAlgError):
            mlab.GaussianKDE(x1, "silverman")

    def test_silverman_singledim_dataset(self):
        """Use a single dimension list as the dataset and test silverman's
        output."""
        x1 = np.array([-7, -5, 1, 4, 5])
        mygauss = mlab.GaussianKDE(x1, "silverman")
        y_expected = 0.76770389927475502
        assert_almost_equal(mygauss.covariance_factor(), y_expected, 7)

    def test_scott_multidim_dataset(self):
        """Use a multi-dimensional array as the dataset and test scott's output
        """
        x1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        with pytest.raises(np.linalg.LinAlgError):
            mlab.GaussianKDE(x1, "scott")

    def test_scott_singledim_dataset(self):
        """Use a single-dimensional array as the dataset and test scott's
        output"""
        x1 = np.array([-7, -5, 1, 4, 5])
        mygauss = mlab.GaussianKDE(x1, "scott")
        y_expected = 0.72477966367769553
        assert_almost_equal(mygauss.covariance_factor(), y_expected, 7)

    def test_scalar_empty_dataset(self):
        """Use an empty array as the dataset and test the scalar's cov factor
        """
        with pytest.raises(ValueError):
            mlab.GaussianKDE([], bw_method=5)

    def test_scalar_covariance_dataset(self):
        """Use a dataset and test a scalar's cov factor
        """
        np.random.seed(8765678)
        n_basesample = 50
        multidim_data = [np.random.randn(n_basesample) for i in range(5)]

        kde = mlab.GaussianKDE(multidim_data, bw_method=0.5)
        assert kde.covariance_factor() == 0.5

    def test_callable_covariance_dataset(self):
        """Use a multi-dimensional array as the dataset and test the callable's
        cov factor"""
        np.random.seed(8765678)
        n_basesample = 50
        multidim_data = [np.random.randn(n_basesample) for i in range(5)]

        def callable_fun(x):
            return 0.55
        kde = mlab.GaussianKDE(multidim_data, bw_method=callable_fun)
        assert kde.covariance_factor() == 0.55

    def test_callable_singledim_dataset(self):
        """Use a single-dimensional array as the dataset and test the
        callable's cov factor"""
        np.random.seed(8765678)
        n_basesample = 50
        multidim_data = np.random.randn(n_basesample)

        kde = mlab.GaussianKDE(multidim_data, bw_method='silverman')
        y_expected = 0.48438841363348911
        assert_almost_equal(kde.covariance_factor(), y_expected, 7)

    def test_wrong_bw_method(self):
        """Test the error message that should be called when bw is invalid."""
        np.random.seed(8765678)
        n_basesample = 50
        data = np.random.randn(n_basesample)
        with pytest.raises(ValueError):
            mlab.GaussianKDE(data, bw_method="invalid")


class TestGaussianKDEEvaluate(object):

    def test_evaluate_diff_dim(self):
        """Test the evaluate method when the dim's of dataset and points are
        different dimensions"""
        x1 = np.arange(3, 10, 2)
        kde = mlab.GaussianKDE(x1)
        x2 = np.arange(3, 12, 2)
        y_expected = [
            0.08797252, 0.11774109, 0.11774109, 0.08797252, 0.0370153
        ]
        y = kde.evaluate(x2)
        np.testing.assert_array_almost_equal(y, y_expected, 7)

    def test_evaluate_inv_dim(self):
        """ Invert the dimensions. i.e., Give the dataset a dimension of
        1 [3,2,4], and the points will have a dimension of 3 [[3],[2],[4]].
        ValueError should be raised"""
        np.random.seed(8765678)
        n_basesample = 50
        multidim_data = np.random.randn(n_basesample)
        kde = mlab.GaussianKDE(multidim_data)
        x2 = [[1], [2], [3]]
        with pytest.raises(ValueError):
            kde.evaluate(x2)

    def test_evaluate_dim_and_num(self):
        """ Tests if evaluated against a one by one array"""
        x1 = np.arange(3, 10, 2)
        x2 = np.array([3])
        kde = mlab.GaussianKDE(x1)
        y_expected = [0.08797252]
        y = kde.evaluate(x2)
        np.testing.assert_array_almost_equal(y, y_expected, 7)

    def test_evaluate_point_dim_not_one(self):
        """Test"""
        x1 = np.arange(3, 10, 2)
        x2 = [np.arange(3, 10, 2), np.arange(3, 10, 2)]
        kde = mlab.GaussianKDE(x1)
        with pytest.raises(ValueError):
            kde.evaluate(x2)

    def test_evaluate_equal_dim_and_num_lt(self):
        """Test when line 3810 fails"""
        x1 = np.arange(3, 10, 2)
        x2 = np.arange(3, 8, 2)
        kde = mlab.GaussianKDE(x1)
        y_expected = [0.08797252, 0.11774109, 0.11774109]
        y = kde.evaluate(x2)
        np.testing.assert_array_almost_equal(y, y_expected, 7)


def test_contiguous_regions():
    a, b, c = 3, 4, 5
    # Starts and ends with True
    mask = [True]*a + [False]*b + [True]*c
    expected = [(0, a), (a+b, a+b+c)]
    with pytest.warns(MatplotlibDeprecationWarning):
        assert mlab.contiguous_regions(mask) == expected
    d, e = 6, 7
    # Starts with True ends with False
    mask = mask + [False]*e
    with pytest.warns(MatplotlibDeprecationWarning):
        assert mlab.contiguous_regions(mask) == expected
    # Starts with False ends with True
    mask = [False]*d + mask[:-e]
    expected = [(d, d+a), (d+a+b, d+a+b+c)]
    with pytest.warns(MatplotlibDeprecationWarning):
        assert mlab.contiguous_regions(mask) == expected
    # Starts and ends with False
    mask = mask + [False]*e
    with pytest.warns(MatplotlibDeprecationWarning):
        assert mlab.contiguous_regions(mask) == expected
        # No True in mask
        assert mlab.contiguous_regions([False]*5) == []
        # Empty mask
        assert mlab.contiguous_regions([]) == []


def test_psd_onesided_norm():
    u = np.array([0, 1, 2, 3, 1, 2, 1])
    dt = 1.0
    Su = np.abs(np.fft.fft(u) * dt)**2 / (dt * u.size)
    P, f = mlab.psd(u, NFFT=u.size, Fs=1/dt, window=mlab.window_none,
                    detrend=mlab.detrend_none, noverlap=0, pad_to=None,
                    scale_by_freq=None,
                    sides='onesided')
    Su_1side = np.append([Su[0]], Su[1:4] + Su[4:][::-1])
    assert_allclose(P, Su_1side, atol=1e-06)


def test_psd_oversampling():
    """Test the case len(x) < NFFT for psd()."""
    u = np.array([0, 1, 2, 3, 1, 2, 1])
    dt = 1.0
    Su = np.abs(np.fft.fft(u) * dt)**2 / (dt * u.size)
    P, f = mlab.psd(u, NFFT=u.size*2, Fs=1/dt, window=mlab.window_none,
                    detrend=mlab.detrend_none, noverlap=0, pad_to=None,
                    scale_by_freq=None,
                    sides='onesided')
    Su_1side = np.append([Su[0]], Su[1:4] + Su[4:][::-1])
    assert_almost_equal(np.sum(P), np.sum(Su_1side))  # same energy
