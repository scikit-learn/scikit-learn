#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

from itertools import combinations
import numpy as np
from numpy.testing import (run_module_suite, assert_almost_equal,
                           assert_allclose, assert_, assert_equal,
                           assert_raises, dec)
import pywt

# Check that float32 and complex64 are preserved.  Other real types get
# converted to float64.
dtypes_in = [np.int8, np.float32, np.float64, np.complex64, np.complex128]
dtypes_out = [np.float64, np.float32, np.float64, np.complex64, np.complex128]


# tolerances used in accuracy comparisons
tol_single = 1e-6
tol_double = 1e-13
dtypes_and_tolerances = [(np.float32, tol_single), (np.float64, tol_double),
                         (np.int8, tol_double), (np.complex64, tol_single),
                         (np.complex128, tol_double)]


# determine which wavelets to test
wavelist = pywt.wavelist()
if 'dmey' in wavelist:
    # accuracy is very low for dmey, so omit it
    wavelist.remove('dmey')
# removing wavelets with dwt_possible == False
del_list = []
for wavelet in wavelist:
    if not isinstance(pywt.DiscreteContinuousWavelet(wavelet), pywt.Wavelet):
        del_list.append(wavelet)
for del_ind in del_list:
    wavelist.remove(del_ind)


####
# 1d multilevel dwt tests
####

def test_wavedec():
    x = [3, 7, 1, 1, -2, 5, 4, 6]
    db1 = pywt.Wavelet('db1')
    cA3, cD3, cD2, cD1 = pywt.wavedec(x, db1)
    assert_almost_equal(cA3, [8.83883476])
    assert_almost_equal(cD3, [-0.35355339])
    assert_allclose(cD2, [4., -3.5])
    assert_allclose(cD1, [-2.82842712, 0, -4.94974747, -1.41421356])
    assert_(pywt.dwt_max_level(len(x), db1) == 3)


def test_waverec_invalid_inputs():
    # input must be list or tuple
    assert_raises(ValueError, pywt.waverec, np.ones(8), 'haar')

    # input list cannot be empty
    assert_raises(ValueError, pywt.waverec, [], 'haar')


def test_waverec_accuracies():
    rstate = np.random.RandomState(1234)
    x0 = rstate.randn(8)
    for dt, tol in dtypes_and_tolerances:
        x = x0.astype(dt)
        if np.iscomplexobj(x):
            x += 1j*rstate.randn(8).astype(x.real.dtype)
        coeffs = pywt.wavedec(x, 'db1')
        assert_allclose(pywt.waverec(coeffs, 'db1'), x, atol=tol, rtol=tol)


def test_waverec_none():
    x = [3, 7, 1, 1, -2, 5, 4, 6]
    coeffs = pywt.wavedec(x, 'db1')

    # set some coefficients to None
    coeffs[2] = None
    coeffs[0] = None
    assert_(pywt.waverec(coeffs, 'db1').size, len(x))


def test_waverec_odd_length():
    x = [3, 7, 1, 1, -2, 5]
    coeffs = pywt.wavedec(x, 'db1')
    assert_allclose(pywt.waverec(coeffs, 'db1'), x, rtol=1e-12)


def test_waverec_complex():
    x = np.array([3, 7, 1, 1, -2, 5, 4, 6])
    x = x + 1j
    coeffs = pywt.wavedec(x, 'db1')
    assert_allclose(pywt.waverec(coeffs, 'db1'), x, rtol=1e-12)


def test_multilevel_dtypes_1d():
    # only checks that the result is of the expected type
    wavelet = pywt.Wavelet('haar')
    for dt_in, dt_out in zip(dtypes_in, dtypes_out):
        # wavedec, waverec
        x = np.ones(8, dtype=dt_in)
        errmsg = "wrong dtype returned for {0} input".format(dt_in)

        coeffs = pywt.wavedec(x, wavelet, level=2)
        for c in coeffs:
            assert_(c.dtype == dt_out, "wavedec: " + errmsg)
        x_roundtrip = pywt.waverec(coeffs, wavelet)
        assert_(x_roundtrip.dtype == dt_out, "waverec: " + errmsg)


def test_waverec_all_wavelets_modes():
    # test 2D case using all wavelets and modes
    rstate = np.random.RandomState(1234)
    r = rstate.randn(80)
    for wavelet in wavelist:
        for mode in pywt.Modes.modes:
            coeffs = pywt.wavedec(r, wavelet, mode=mode)
            assert_allclose(pywt.waverec(coeffs, wavelet, mode=mode),
                            r, rtol=tol_single, atol=tol_single)

####
# 2d multilevel dwt function tests
####


def test_waverec2_accuracies():
    rstate = np.random.RandomState(1234)
    x0 = rstate.randn(4, 4)
    for dt, tol in dtypes_and_tolerances:
        x = x0.astype(dt)
        if np.iscomplexobj(x):
            x += 1j*rstate.randn(4, 4).astype(x.real.dtype)
        coeffs = pywt.wavedec2(x, 'db1')
        assert_(len(coeffs) == 3)
        assert_allclose(pywt.waverec2(coeffs, 'db1'), x, atol=tol, rtol=tol)


def test_multilevel_dtypes_2d():
    wavelet = pywt.Wavelet('haar')
    for dt_in, dt_out in zip(dtypes_in, dtypes_out):
        # wavedec2, waverec2
        x = np.ones((8, 8), dtype=dt_in)
        errmsg = "wrong dtype returned for {0} input".format(dt_in)
        cA, coeffsD2, coeffsD1 = pywt.wavedec2(x, wavelet, level=2)
        assert_(cA.dtype == dt_out, "wavedec2: " + errmsg)
        for c in coeffsD1:
            assert_(c.dtype == dt_out, "wavedec2: " + errmsg)
        for c in coeffsD2:
            assert_(c.dtype == dt_out, "wavedec2: " + errmsg)
        x_roundtrip = pywt.waverec2([cA, coeffsD2, coeffsD1], wavelet)
        assert_(x_roundtrip.dtype == dt_out, "waverec2: " + errmsg)


@dec.slow
def test_waverec2_all_wavelets_modes():
    # test 2D case using all wavelets and modes
    rstate = np.random.RandomState(1234)
    r = rstate.randn(80, 96)
    for wavelet in wavelist:
        for mode in pywt.Modes.modes:
            coeffs = pywt.wavedec2(r, wavelet, mode=mode)
            assert_allclose(pywt.waverec2(coeffs, wavelet, mode=mode),
                            r, rtol=tol_single, atol=tol_single)


def test_wavedec2_complex():
    data = np.ones((4, 4)) + 1j
    coeffs = pywt.wavedec2(data, 'db1')
    assert_(len(coeffs) == 3)
    assert_allclose(pywt.waverec2(coeffs, 'db1'), data, rtol=1e-12)


def test_wavedec2_invalid_inputs():
    # input array has too few dimensions
    data = np.ones(4)
    assert_raises(ValueError, pywt.wavedec2, data, 'haar')


def test_waverec2_invalid_inputs():
    # input must be list or tuple
    assert_raises(ValueError, pywt.waverec2, np.ones((8, 8)), 'haar')

    # input list cannot be empty
    assert_raises(ValueError, pywt.waverec2, [], 'haar')


def test_waverec2_coeff_shape_mismatch():
    x = np.ones((8, 8))
    coeffs = pywt.wavedec2(x, 'db1')

    # introduce a shape mismatch in the coefficients
    coeffs = list(coeffs)
    coeffs[1] = list(coeffs[1])
    coeffs[1][1] = np.zeros((16, 1))
    assert_raises(ValueError, pywt.waverec2, coeffs, 'db1')


def test_waverec2_odd_length():
    x = np.ones((10, 6))
    coeffs = pywt.wavedec2(x, 'db1')
    assert_allclose(pywt.waverec2(coeffs, 'db1'), x, rtol=1e-12)


def test_waverec2_none_coeffs():
    x = np.arange(24).reshape(6, 4)
    coeffs = pywt.wavedec2(x, 'db1')
    coeffs[1] = (None, None, None)
    assert_(x.shape == pywt.waverec2(coeffs, 'db1').shape)

####
# nd multilevel dwt function tests
####


def test_waverecn():
    rstate = np.random.RandomState(1234)
    # test 1D through 4D cases
    for nd in range(1, 5):
        x = rstate.randn(*(4, )*nd)
        coeffs = pywt.wavedecn(x, 'db1')
        assert_(len(coeffs) == 3)
        assert_allclose(pywt.waverecn(coeffs, 'db1'), x, rtol=tol_double)


def test_waverecn_empty_coeff():
    coeffs = [np.ones((2, 2, 2)), {}, {}]
    assert_equal(pywt.waverecn(coeffs, 'db1').shape, (8, 8, 8))

    assert_equal(pywt.waverecn(coeffs, 'db1').shape, (8, 8, 8))
    coeffs = [np.ones((2, 2, 2)), {}, {'daa': np.ones((4, 4, 4))}]

    coeffs = [np.ones((2, 2, 2)), {}, {}, {'daa': np.ones((8, 8, 8))}]
    assert_equal(pywt.waverecn(coeffs, 'db1').shape, (16, 16, 16))


def test_waverecn_invalid_coeffs():
    # approximation coeffs as None and no valid detail oeffs
    coeffs = [None, {}]
    assert_raises(ValueError, pywt.waverecn, coeffs, 'db1')

    # use of None for a coefficient value
    coeffs = [np.ones((2, 2, 2)), {}, {'daa': None}, ]
    assert_raises(ValueError, pywt.waverecn, coeffs, 'db1')

    # invalid key names in coefficient list
    coeffs = [np.ones((4, 4, 4)), {'daa': np.ones((4, 4, 4)),
                                   'foo': np.ones((4, 4, 4))}]
    assert_raises(ValueError, pywt.waverecn, coeffs, 'db1')

    # mismatched key name lengths
    coeffs = [np.ones((4, 4, 4)), {'daa': np.ones((4, 4, 4)),
                                   'da': np.ones((4, 4, 4))}]
    assert_raises(ValueError, pywt.waverecn, coeffs, 'db1')

    # key name lengths don't match the array dimensions
    coeffs = [[[[1.0]]], {'ad': [[[0.0]]], 'da': [[[0.0]]], 'dd': [[[0.0]]]}]
    assert_raises(ValueError, pywt.waverecn, coeffs, 'db1')

    # input list cannot be empty
    assert_raises(ValueError, pywt.waverecn, [], 'haar')


def test_waverecn_lists():
    # support coefficient arrays specified as lists instead of arrays
    coeffs = [[[1.0]], {'ad': [[0.0]], 'da': [[0.0]], 'dd': [[0.0]]}]
    assert_equal(pywt.waverecn(coeffs, 'db1').shape, (2, 2))


def test_waverecn_invalid_coeffs2():
    # shape mismatch should raise an error
    coeffs = [np.ones((4, 4, 4)), {'ada': np.ones((4, 4))}]
    assert_raises(ValueError, pywt.waverecn, coeffs, 'db1')


def test_wavedecn_invalid_inputs():
    # input array has too few dimensions
    data = np.array(0)
    assert_raises(ValueError, pywt.wavedecn, data, 'haar')

    # invalid number of levels
    data = np.ones(16)
    assert_raises(ValueError, pywt.wavedecn, data, 'haar', level=-1)
    assert_raises(ValueError, pywt.wavedecn, data, 'haar', level=100)


def test_waverecn_accuracies():
    # testing 3D only here
    rstate = np.random.RandomState(1234)
    x0 = rstate.randn(4, 4, 4)
    for dt, tol in dtypes_and_tolerances:
        x = x0.astype(dt)
        if np.iscomplexobj(x):
            x += 1j*rstate.randn(4, 4, 4).astype(x.real.dtype)
        coeffs = pywt.wavedecn(x.astype(dt), 'db1')
        assert_allclose(pywt.waverecn(coeffs, 'db1'), x, atol=tol, rtol=tol)


def test_multilevel_dtypes_nd():
    wavelet = pywt.Wavelet('haar')
    for dt_in, dt_out in zip(dtypes_in, dtypes_out):
        # wavedecn, waverecn
        x = np.ones((8, 8), dtype=dt_in)
        errmsg = "wrong dtype returned for {0} input".format(dt_in)
        cA, coeffsD2, coeffsD1 = pywt.wavedecn(x, wavelet, level=2)
        assert_(cA.dtype == dt_out, "wavedecn: " + errmsg)
        for key, c in coeffsD1.items():
            assert_(c.dtype == dt_out, "wavedecn: " + errmsg)
        for key, c in coeffsD2.items():
            assert_(c.dtype == dt_out, "wavedecn: " + errmsg)
        x_roundtrip = pywt.waverecn([cA, coeffsD2, coeffsD1], wavelet)
        assert_(x_roundtrip.dtype == dt_out, "waverecn: " + errmsg)


def test_wavedecn_complex():
    data = np.ones((4, 4, 4)) + 1j
    coeffs = pywt.wavedecn(data, 'db1')
    assert_allclose(pywt.waverecn(coeffs, 'db1'), data, rtol=1e-12)


def test_waverecn_dtypes():
    x = np.ones((4, 4, 4))
    for dt, tol in dtypes_and_tolerances:
        coeffs = pywt.wavedecn(x.astype(dt), 'db1')
        assert_allclose(pywt.waverecn(coeffs, 'db1'), x, atol=tol, rtol=tol)


@dec.slow
def test_waverecn_all_wavelets_modes():
    # test 2D case using all wavelets and modes
    rstate = np.random.RandomState(1234)
    r = rstate.randn(80, 96)
    for wavelet in wavelist:
        for mode in pywt.Modes.modes:
            coeffs = pywt.wavedecn(r, wavelet, mode=mode)
            assert_allclose(pywt.waverecn(coeffs, wavelet, mode=mode),
                            r, rtol=tol_single, atol=tol_single)


def test_coeffs_to_array():
    # single element list returns the first element
    a_coeffs = [np.arange(8).reshape(2, 4), ]
    arr, arr_slices = pywt.coeffs_to_array(a_coeffs)
    assert_allclose(arr, a_coeffs[0])
    assert_allclose(arr, arr[arr_slices[0]])

    assert_raises(ValueError, pywt.coeffs_to_array, [])
    # invalid second element:  array as in wavedec, but not 1D
    assert_raises(ValueError, pywt.coeffs_to_array, [a_coeffs[0], ] * 2)
    # invalid second element:  tuple as in wavedec2, but not a 3-tuple
    assert_raises(ValueError, pywt.coeffs_to_array, [a_coeffs[0],
                                                     (a_coeffs[0], )])
    # coefficients as None is not supported
    assert_raises(ValueError, pywt.coeffs_to_array, [None, ])
    assert_raises(ValueError, pywt.coeffs_to_array, [a_coeffs,
                                                     (None, None, None)])

    # invalid type for second coefficient list element
    assert_raises(ValueError, pywt.coeffs_to_array, [a_coeffs, None])

    # use an invalid key name in the coef dictionary
    coeffs = [np.array([0]), dict(d=np.array([0]), c=np.array([0]))]
    assert_raises(ValueError, pywt.coeffs_to_array, coeffs)


def test_wavedecn_coeff_reshape_even():
    # verify round trip is correct:
    #   wavedecn - >coeffs_to_array-> array_to_coeffs -> waverecn
    # This is done for wavedec{1, 2, n}
    rng = np.random.RandomState(1234)
    params = {'wavedec': {'d': 1, 'dec': pywt.wavedec, 'rec': pywt.waverec},
              'wavedec2': {'d': 2, 'dec': pywt.wavedec2, 'rec': pywt.waverec2},
              'wavedecn': {'d': 3, 'dec': pywt.wavedecn, 'rec': pywt.waverecn}}
    N = 28
    for f in params:
        x1 = rng.randn(*([N] * params[f]['d']))
        for mode in pywt.Modes.modes:
            for wave in wavelist:
                w = pywt.Wavelet(wave)
                maxlevel = pywt.dwt_max_level(np.min(x1.shape), w.dec_len)
                if maxlevel == 0:
                    continue

                coeffs = params[f]['dec'](x1, w, mode=mode)
                coeff_arr, coeff_slices = pywt.coeffs_to_array(coeffs)
                coeffs2 = pywt.array_to_coeffs(coeff_arr, coeff_slices,
                                               output_format=f)
                x1r = params[f]['rec'](coeffs2, w, mode=mode)

                assert_allclose(x1, x1r, rtol=1e-4, atol=1e-4)


def test_wavedecn_coeff_reshape_axes_subset():
    # verify round trip is correct when only a subset of axes are transformed:
    #   wavedecn - >coeffs_to_array-> array_to_coeffs -> waverecn
    # This is done for wavedec{1, 2, n}
    rng = np.random.RandomState(1234)
    mode = 'symmetric'
    w = pywt.Wavelet('db2')
    N = 16
    ndim = 3
    for axes in [(-1, ), (0, ), (1, ), (0, 1), (1, 2), (0, 2), None]:
        x1 = rng.randn(*([N] * ndim))
        coeffs = pywt.wavedecn(x1, w, mode=mode, axes=axes)
        coeff_arr, coeff_slices = pywt.coeffs_to_array(coeffs, axes=axes)
        if axes is not None:
            # if axes is not None, it must be provided to coeffs_to_array
            assert_raises(ValueError, pywt.coeffs_to_array, coeffs)

        # mismatched axes size
        assert_raises(ValueError, pywt.coeffs_to_array, coeffs,
                      axes=(0, 1, 2, 3))
        assert_raises(ValueError, pywt.coeffs_to_array, coeffs,
                      axes=())

        coeffs2 = pywt.array_to_coeffs(coeff_arr, coeff_slices)
        x1r = pywt.waverecn(coeffs2, w, mode=mode, axes=axes)

        assert_allclose(x1, x1r, rtol=1e-4, atol=1e-4)


def test_coeffs_to_array_padding():
    rng = np.random.RandomState(1234)
    x1 = rng.randn(32, 32)
    mode = 'symmetric'
    coeffs = pywt.wavedecn(x1, 'db2', mode=mode)

    # padding=None raises a ValueError when tight packing is not possible
    assert_raises(ValueError, pywt.coeffs_to_array, coeffs, padding=None)

    # set padded values to nan
    coeff_arr, coeff_slices = pywt.coeffs_to_array(coeffs, padding=np.nan)
    npad = np.sum(np.isnan(coeff_arr))
    assert_(npad > 0)

    # pad with zeros
    coeff_arr, coeff_slices = pywt.coeffs_to_array(coeffs, padding=0)
    assert_(np.sum(np.isnan(coeff_arr)) == 0)
    assert_(np.sum(coeff_arr == 0) == npad)

    # Haar case with N as a power of 2 can be tightly packed
    coeffs_haar = pywt.wavedecn(x1, 'haar', mode=mode)
    coeff_arr, coeff_slices = pywt.coeffs_to_array(coeffs_haar, padding=None)
    # shape of coeff_arr will match in this case, but not in general
    assert_equal(coeff_arr.shape, x1.shape)


def test_waverecn_coeff_reshape_odd():
    # verify round trip is correct:
    #   wavedecn - >coeffs_to_array-> array_to_coeffs -> waverecn
    rng = np.random.RandomState(1234)
    x1 = rng.randn(35, 33)
    for mode in pywt.Modes.modes:
        for wave in ['haar', ]:
            w = pywt.Wavelet(wave)
            maxlevel = pywt.dwt_max_level(np.min(x1.shape), w.dec_len)
            if maxlevel == 0:
                continue
            coeffs = pywt.wavedecn(x1, w, mode=mode)
            coeff_arr, coeff_slices = pywt.coeffs_to_array(coeffs)
            coeffs2 = pywt.array_to_coeffs(coeff_arr, coeff_slices)
            x1r = pywt.waverecn(coeffs2, w, mode=mode)
            # truncate reconstructed values to original shape
            x1r = x1r[[slice(s) for s in x1.shape]]
            assert_allclose(x1, x1r, rtol=1e-4, atol=1e-4)


def test_array_to_coeffs_invalid_inputs():
    coeffs = pywt.wavedecn(np.ones(2), 'haar')
    arr, arr_slices = pywt.coeffs_to_array(coeffs)

    # empty list of array slices
    assert_raises(ValueError, pywt.array_to_coeffs, arr, [])

    # invalid format name
    assert_raises(ValueError, pywt.array_to_coeffs, arr, arr_slices, 'foo')


def test_waverec_axes_subsets():
    rstate = np.random.RandomState(0)
    data = rstate.standard_normal((8, 8, 8))
    # test all combinations of 1 out of 3 axes transformed
    for axis in [0, 1, 2]:
        coefs = pywt.wavedec(data, 'haar', axis=axis)
        rec = pywt.waverec(coefs, 'haar', axis=axis)
        assert_allclose(rec, data, atol=1e-14)


def test_waverec_axis_db2():
    # test for fix to issue gh-293
    rstate = np.random.RandomState(0)
    data = rstate.standard_normal((16, 16))
    for axis in [0, 1]:
        coefs = pywt.wavedec(data, 'db2', axis=axis)
        rec = pywt.waverec(coefs, 'db2', axis=axis)
        assert_allclose(rec, data, atol=1e-14)


def test_waverec2_axes_subsets():
    rstate = np.random.RandomState(0)
    data = rstate.standard_normal((8, 8, 8))
    # test all combinations of 2 out of 3 axes transformed
    for axes in combinations((0, 1, 2), 2):
        coefs = pywt.wavedec2(data, 'haar', axes=axes)
        rec = pywt.waverec2(coefs, 'haar', axes=axes)
        assert_allclose(rec, data, atol=1e-14)


def test_waverecn_axes_subsets():
    rstate = np.random.RandomState(0)
    data = rstate.standard_normal((8, 8, 8, 8))
    # test all combinations of 3 out of 4 axes transformed
    for axes in combinations((0, 1, 2, 3), 3):
        coefs = pywt.wavedecn(data, 'haar', axes=axes)
        rec = pywt.waverecn(coefs, 'haar', axes=axes)
        assert_allclose(rec, data, atol=1e-14)


def test_waverecn_int_axis():
    # waverecn should also work for axes as an integer
    rstate = np.random.RandomState(0)
    data = rstate.standard_normal((8, 8))
    for axis in [0, 1]:
        coefs = pywt.wavedecn(data, 'haar', axes=axis)
        rec = pywt.waverecn(coefs, 'haar', axes=axis)
        assert_allclose(rec, data, atol=1e-14)


def test_wavedec_axis_error():
    data = np.ones(4)
    # out of range axis not allowed
    assert_raises(ValueError, pywt.wavedec, data, 'haar', axis=1)


def test_waverec_axis_error():
    c = pywt.wavedec(np.ones(4), 'haar')
    # out of range axis not allowed
    assert_raises(ValueError, pywt.waverec, c, 'haar', axis=1)


def test_waverec_shape_mismatch_error():
    c = pywt.wavedec(np.ones(16), 'haar')
    # truncate a detail coefficient to an incorrect shape
    c[3] = c[3][:-1]
    assert_raises(ValueError, pywt.waverec, c, 'haar', axis=1)


def test_wavedec2_axes_errors():
    data = np.ones((4, 4))
    # integer axes not allowed
    assert_raises(TypeError, pywt.wavedec2, data, 'haar', axes=1)
    # non-unique axes not allowed
    assert_raises(ValueError, pywt.wavedec2, data, 'haar', axes=(0, 0))
    # out of range axis not allowed
    assert_raises(ValueError, pywt.wavedec2, data, 'haar', axes=(0, 2))


def test_waverec2_axes_errors():
    data = np.ones((4, 4))
    c = pywt.wavedec2(data, 'haar')
    # integer axes not allowed
    assert_raises(TypeError, pywt.waverec2, c, 'haar', axes=1)
    # non-unique axes not allowed
    assert_raises(ValueError, pywt.waverec2, c, 'haar', axes=(0, 0))
    # out of range axis not allowed
    assert_raises(ValueError, pywt.waverec2, c, 'haar', axes=(0, 2))


def test_wavedecn_axes_errors():
    data = np.ones((8, 8, 8))
    # repeated axes not allowed
    assert_raises(ValueError, pywt.wavedecn, data, 'haar', axes=(1, 1))
    # out of range axis not allowed
    assert_raises(ValueError, pywt.wavedecn, data, 'haar', axes=(0, 1, 3))


def test_waverecn_axes_errors():
    data = np.ones((8, 8, 8))
    c = pywt.wavedecn(data, 'haar')
    # repeated axes not allowed
    assert_raises(ValueError, pywt.waverecn, c, 'haar', axes=(1, 1))
    # out of range axis not allowed
    assert_raises(ValueError, pywt.waverecn, c, 'haar', axes=(0, 1, 3))


if __name__ == '__main__':
    run_module_suite()
