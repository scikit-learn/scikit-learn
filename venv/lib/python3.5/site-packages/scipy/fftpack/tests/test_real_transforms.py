from __future__ import division, print_function, absolute_import

from os.path import join, dirname

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_equal
from pytest import raises as assert_raises

from scipy.fftpack.realtransforms import (
    dct, idct, dst, idst, dctn, idctn, dstn, idstn)

# Matlab reference data
MDATA = np.load(join(dirname(__file__), 'test.npz'))
X = [MDATA['x%d' % i] for i in range(8)]
Y = [MDATA['y%d' % i] for i in range(8)]

# FFTW reference data: the data are organized as follows:
#    * SIZES is an array containing all available sizes
#    * for every type (1, 2, 3, 4) and every size, the array dct_type_size
#    contains the output of the DCT applied to the input np.linspace(0, size-1,
#    size)
FFTWDATA_DOUBLE = np.load(join(dirname(__file__), 'fftw_double_ref.npz'))
FFTWDATA_SINGLE = np.load(join(dirname(__file__), 'fftw_single_ref.npz'))
FFTWDATA_SIZES = FFTWDATA_DOUBLE['sizes']


def fftw_dct_ref(type, size, dt):
    x = np.linspace(0, size-1, size).astype(dt)
    dt = np.result_type(np.float32, dt)
    if dt == np.double:
        data = FFTWDATA_DOUBLE
    elif dt == np.float32:
        data = FFTWDATA_SINGLE
    else:
        raise ValueError()
    y = (data['dct_%d_%d' % (type, size)]).astype(dt)
    return x, y, dt


def fftw_dst_ref(type, size, dt):
    x = np.linspace(0, size-1, size).astype(dt)
    dt = np.result_type(np.float32, dt)
    if dt == np.double:
        data = FFTWDATA_DOUBLE
    elif dt == np.float32:
        data = FFTWDATA_SINGLE
    else:
        raise ValueError()
    y = (data['dst_%d_%d' % (type, size)]).astype(dt)
    return x, y, dt


def dct_2d_ref(x, **kwargs):
    """ used as a reference in testing dct2. """
    x = np.array(x, copy=True)
    for row in range(x.shape[0]):
        x[row, :] = dct(x[row, :], **kwargs)
    for col in range(x.shape[1]):
        x[:, col] = dct(x[:, col], **kwargs)
    return x


def idct_2d_ref(x, **kwargs):
    """ used as a reference in testing idct2. """
    x = np.array(x, copy=True)
    for row in range(x.shape[0]):
        x[row, :] = idct(x[row, :], **kwargs)
    for col in range(x.shape[1]):
        x[:, col] = idct(x[:, col], **kwargs)
    return x


def dst_2d_ref(x, **kwargs):
    """ used as a reference in testing dst2. """
    x = np.array(x, copy=True)
    for row in range(x.shape[0]):
        x[row, :] = dst(x[row, :], **kwargs)
    for col in range(x.shape[1]):
        x[:, col] = dst(x[:, col], **kwargs)
    return x


def idst_2d_ref(x, **kwargs):
    """ used as a reference in testing idst2. """
    x = np.array(x, copy=True)
    for row in range(x.shape[0]):
        x[row, :] = idst(x[row, :], **kwargs)
    for col in range(x.shape[1]):
        x[:, col] = idst(x[:, col], **kwargs)
    return x


class TestComplex(object):
    def test_dct_complex64(self):
        y = dct(1j*np.arange(5, dtype=np.complex64))
        x = 1j*dct(np.arange(5))
        assert_array_almost_equal(x, y)

    def test_dct_complex(self):
        y = dct(np.arange(5)*1j)
        x = 1j*dct(np.arange(5))
        assert_array_almost_equal(x, y)

    def test_idct_complex(self):
        y = idct(np.arange(5)*1j)
        x = 1j*idct(np.arange(5))
        assert_array_almost_equal(x, y)

    def test_dst_complex64(self):
        y = dst(np.arange(5, dtype=np.complex64)*1j)
        x = 1j*dst(np.arange(5))
        assert_array_almost_equal(x, y)

    def test_dst_complex(self):
        y = dst(np.arange(5)*1j)
        x = 1j*dst(np.arange(5))
        assert_array_almost_equal(x, y)

    def test_idst_complex(self):
        y = idst(np.arange(5)*1j)
        x = 1j*idst(np.arange(5))
        assert_array_almost_equal(x, y)


class _TestDCTBase(object):
    def setup_method(self):
        self.rdt = None
        self.dec = 14
        self.type = None

    def test_definition(self):
        for i in FFTWDATA_SIZES:
            x, yr, dt = fftw_dct_ref(self.type, i, self.rdt)
            y = dct(x, type=self.type)
            assert_equal(y.dtype, dt)
            # XXX: we divide by np.max(y) because the tests fail otherwise. We
            # should really use something like assert_array_approx_equal. The
            # difference is due to fftw using a better algorithm w.r.t error
            # propagation compared to the ones from fftpack.
            assert_array_almost_equal(y / np.max(y), yr / np.max(y), decimal=self.dec,
                    err_msg="Size %d failed" % i)

    def test_axis(self):
        nt = 2
        for i in [7, 8, 9, 16, 32, 64]:
            x = np.random.randn(nt, i)
            y = dct(x, type=self.type)
            for j in range(nt):
                assert_array_almost_equal(y[j], dct(x[j], type=self.type),
                        decimal=self.dec)

            x = x.T
            y = dct(x, axis=0, type=self.type)
            for j in range(nt):
                assert_array_almost_equal(y[:,j], dct(x[:,j], type=self.type),
                        decimal=self.dec)


class _TestDCTIIBase(_TestDCTBase):
    def test_definition_matlab(self):
        # Test correspondence with matlab (orthornomal mode).
        for i in range(len(X)):
            dt = np.result_type(np.float32, self.rdt)
            x = np.array(X[i], dtype=dt)

            yr = Y[i]
            y = dct(x, norm="ortho", type=2)
            assert_equal(y.dtype, dt)
            assert_array_almost_equal(y, yr, decimal=self.dec)


class _TestDCTIIIBase(_TestDCTBase):
    def test_definition_ortho(self):
        # Test orthornomal mode.
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            dt = np.result_type(np.float32, self.rdt)
            y = dct(x, norm='ortho', type=2)
            xi = dct(y, norm="ortho", type=3)
            assert_equal(xi.dtype, dt)
            assert_array_almost_equal(xi, x, decimal=self.dec)


class TestDCTIDouble(_TestDCTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 1


class TestDCTIFloat(_TestDCTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 1


class TestDCTIInt(_TestDCTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 5
        self.type = 1


class TestDCTIIDouble(_TestDCTIIBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 2


class TestDCTIIFloat(_TestDCTIIBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 2


class TestDCTIIInt(_TestDCTIIBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 5
        self.type = 2


class TestDCTIIIDouble(_TestDCTIIIBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 3


class TestDCTIIIFloat(_TestDCTIIIBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 3


class TestDCTIIIInt(_TestDCTIIIBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 5
        self.type = 3


class _TestIDCTBase(object):
    def setup_method(self):
        self.rdt = None
        self.dec = 14
        self.type = None

    def test_definition(self):
        for i in FFTWDATA_SIZES:
            xr, yr, dt = fftw_dct_ref(self.type, i, self.rdt)
            x = idct(yr, type=self.type)
            if self.type == 1:
                x /= 2 * (i-1)
            else:
                x /= 2 * i
            assert_equal(x.dtype, dt)
            # XXX: we divide by np.max(y) because the tests fail otherwise. We
            # should really use something like assert_array_approx_equal. The
            # difference is due to fftw using a better algorithm w.r.t error
            # propagation compared to the ones from fftpack.
            assert_array_almost_equal(x / np.max(x), xr / np.max(x), decimal=self.dec,
                    err_msg="Size %d failed" % i)


class TestIDCTIDouble(_TestIDCTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 1


class TestIDCTIFloat(_TestIDCTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 4
        self.type = 1


class TestIDCTIInt(_TestIDCTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 4
        self.type = 1


class TestIDCTIIDouble(_TestIDCTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 2


class TestIDCTIIFloat(_TestIDCTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 2


class TestIDCTIIInt(_TestIDCTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 5
        self.type = 2


class TestIDCTIIIDouble(_TestIDCTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 3


class TestIDCTIIIFloat(_TestIDCTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 3


class TestIDCTIIIInt(_TestIDCTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 5
        self.type = 3


class _TestDSTBase(object):
    def setup_method(self):
        self.rdt = None  # dtype
        self.dec = None  # number of decimals to match
        self.type = None  # dst type

    def test_definition(self):
        for i in FFTWDATA_SIZES:
            xr, yr, dt = fftw_dst_ref(self.type, i, self.rdt)
            y = dst(xr, type=self.type)
            assert_equal(y.dtype, dt)
            # XXX: we divide by np.max(y) because the tests fail otherwise. We
            # should really use something like assert_array_approx_equal. The
            # difference is due to fftw using a better algorithm w.r.t error
            # propagation compared to the ones from fftpack.
            assert_array_almost_equal(y / np.max(y), yr / np.max(y), decimal=self.dec,
                    err_msg="Size %d failed" % i)


class TestDSTIDouble(_TestDSTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 1


class TestDSTIFloat(_TestDSTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 1


class TestDSTIInt(_TestDSTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 5
        self.type = 1


class TestDSTIIDouble(_TestDSTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 2


class TestDSTIIFloat(_TestDSTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 6
        self.type = 2


class TestDSTIIInt(_TestDSTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 6
        self.type = 2


class TestDSTIIIDouble(_TestDSTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 3


class TestDSTIIIFloat(_TestDSTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 7
        self.type = 3


class TestDSTIIIInt(_TestDSTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 7
        self.type = 3


class _TestIDSTBase(object):
    def setup_method(self):
        self.rdt = None
        self.dec = None
        self.type = None

    def test_definition(self):
        for i in FFTWDATA_SIZES:
            xr, yr, dt = fftw_dst_ref(self.type, i, self.rdt)
            x = idst(yr, type=self.type)
            if self.type == 1:
                x /= 2 * (i+1)
            else:
                x /= 2 * i
            assert_equal(x.dtype, dt)
            # XXX: we divide by np.max(x) because the tests fail otherwise. We
            # should really use something like assert_array_approx_equal. The
            # difference is due to fftw using a better algorithm w.r.t error
            # propagation compared to the ones from fftpack.
            assert_array_almost_equal(x / np.max(x), xr / np.max(x), decimal=self.dec,
                    err_msg="Size %d failed" % i)


class TestIDSTIDouble(_TestIDSTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 12
        self.type = 1


class TestIDSTIFloat(_TestIDSTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 4
        self.type = 1


class TestIDSTIInt(_TestIDSTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 4
        self.type = 1


class TestIDSTIIDouble(_TestIDSTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 2


class TestIDSTIIFloat(_TestIDSTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 6
        self.type = 2


class TestIDSTIIInt(_TestIDSTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 6
        self.type = 2


class TestIDSTIIIDouble(_TestIDSTBase):
    def setup_method(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 3


class TestIDSTIIIFloat(_TestIDSTBase):
    def setup_method(self):
        self.rdt = np.float32
        self.dec = 6
        self.type = 3


class TestIDSTIIIInt(_TestIDSTBase):
    def setup_method(self):
        self.rdt = int
        self.dec = 6
        self.type = 3


class TestOverwrite(object):
    """Check input overwrite behavior """

    real_dtypes = [np.float32, np.float64]

    def _check(self, x, routine, type, fftsize, axis, norm, overwrite_x,
               should_overwrite, **kw):
        x2 = x.copy()
        routine(x2, type, fftsize, axis, norm, overwrite_x=overwrite_x)

        sig = "%s(%s%r, %r, axis=%r, overwrite_x=%r)" % (
            routine.__name__, x.dtype, x.shape, fftsize, axis, overwrite_x)
        if not should_overwrite:
            assert_equal(x2, x, err_msg="spurious overwrite in %s" % sig)

    def _check_1d(self, routine, dtype, shape, axis, overwritable_dtypes):
        np.random.seed(1234)
        if np.issubdtype(dtype, np.complexfloating):
            data = np.random.randn(*shape) + 1j*np.random.randn(*shape)
        else:
            data = np.random.randn(*shape)
        data = data.astype(dtype)

        for type in [1, 2, 3]:
            for overwrite_x in [True, False]:
                for norm in [None, 'ortho']:
                    if type == 1 and norm == 'ortho':
                        continue

                    should_overwrite = (overwrite_x
                                        and dtype in overwritable_dtypes
                                        and (len(shape) == 1 or
                                             (axis % len(shape) == len(shape)-1
                                              )))
                    self._check(data, routine, type, None, axis, norm,
                                overwrite_x, should_overwrite)

    def test_dct(self):
        overwritable = self.real_dtypes
        for dtype in self.real_dtypes:
            self._check_1d(dct, dtype, (16,), -1, overwritable)
            self._check_1d(dct, dtype, (16, 2), 0, overwritable)
            self._check_1d(dct, dtype, (2, 16), 1, overwritable)

    def test_idct(self):
        overwritable = self.real_dtypes
        for dtype in self.real_dtypes:
            self._check_1d(idct, dtype, (16,), -1, overwritable)
            self._check_1d(idct, dtype, (16, 2), 0, overwritable)
            self._check_1d(idct, dtype, (2, 16), 1, overwritable)

    def test_dst(self):
        overwritable = self.real_dtypes
        for dtype in self.real_dtypes:
            self._check_1d(dst, dtype, (16,), -1, overwritable)
            self._check_1d(dst, dtype, (16, 2), 0, overwritable)
            self._check_1d(dst, dtype, (2, 16), 1, overwritable)

    def test_idst(self):
        overwritable = self.real_dtypes
        for dtype in self.real_dtypes:
            self._check_1d(idst, dtype, (16,), -1, overwritable)
            self._check_1d(idst, dtype, (16, 2), 0, overwritable)
            self._check_1d(idst, dtype, (2, 16), 1, overwritable)


class Test_DCTN_IDCTN(object):
    dec = 14
    types = [1, 2, 3]
    norms = [None, 'ortho']
    rstate = np.random.RandomState(1234)
    shape = (32, 16)
    data = rstate.randn(*shape)
    # Sets of functions to test
    function_sets = [dict(forward=dctn,
                          inverse=idctn,
                          forward_ref=dct_2d_ref,
                          inverse_ref=idct_2d_ref),
                     dict(forward=dstn,
                          inverse=idstn,
                          forward_ref=dst_2d_ref,
                          inverse_ref=idst_2d_ref), ]

    def test_axes_round_trip(self):
        norm = 'ortho'
        for function_set in self.function_sets:
            fforward = function_set['forward']
            finverse = function_set['inverse']
            for axes in [None, (1, ), (0, ), (0, 1), (-2, -1)]:
                for dct_type in self.types:
                    if norm == 'ortho' and dct_type == 1:
                        continue  # 'ortho' not supported by DCT-I
                    tmp = fforward(self.data, type=dct_type, axes=axes,
                                   norm=norm)
                    tmp = finverse(tmp, type=dct_type, axes=axes, norm=norm)
                    assert_array_almost_equal(self.data, tmp, decimal=self.dec)

    def test_dctn_vs_2d_reference(self):
        for function_set in self.function_sets:
            fforward = function_set['forward']
            fforward_ref = function_set['forward_ref']
            for dct_type in self.types:
                for norm in self.norms:
                    if norm == 'ortho' and dct_type == 1:
                        continue  # 'ortho' not supported by DCT-I
                    y1 = fforward(self.data, type=dct_type, axes=None,
                                  norm=norm)
                    y2 = fforward_ref(self.data, type=dct_type, norm=norm)
                    assert_array_almost_equal(y1, y2, decimal=11)

    def test_idctn_vs_2d_reference(self):
        for function_set in self.function_sets:
            finverse = function_set['inverse']
            finverse_ref = function_set['inverse_ref']
            for dct_type in self.types:
                for norm in self.norms:
                    print(function_set, dct_type, norm)
                    if norm == 'ortho' and dct_type == 1:
                        continue  # 'ortho' not supported by DCT-I
                    fdata = dctn(self.data, type=dct_type, norm=norm)
                    y1 = finverse(fdata, type=dct_type, norm=norm)
                    y2 = finverse_ref(fdata, type=dct_type, norm=norm)
                    assert_array_almost_equal(y1, y2, decimal=11)

    def test_axes_and_shape(self):
        for function_set in self.function_sets:
            fforward = function_set['forward']
            finverse = function_set['inverse']

            # shape must match the number of axes
            assert_raises(ValueError, fforward, self.data,
                          shape=(self.data.shape[0], ),
                          axes=(0, 1))
            assert_raises(ValueError, fforward, self.data,
                          shape=(self.data.shape[0], ),
                          axes=None)
            assert_raises(ValueError, fforward, self.data,
                          shape=self.data.shape,
                          axes=(0, ))
            # shape must be a tuple
            assert_raises(TypeError, fforward, self.data,
                          shape=self.data.shape[0],
                          axes=(0, 1))

            # shape=None works with a subset of axes
            for axes in [(0, ), (1, )]:
                tmp = fforward(self.data, shape=None, axes=axes, norm='ortho')
                tmp = finverse(tmp, shape=None, axes=axes, norm='ortho')
                assert_array_almost_equal(self.data, tmp, decimal=self.dec)

            # non-default shape
            tmp = fforward(self.data, shape=(128, 128), axes=None)
            assert_equal(tmp.shape, (128, 128))
