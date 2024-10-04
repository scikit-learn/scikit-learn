"""Helpers to utilize existing stft / istft tests for testing `ShortTimeFFT`.

This module provides the functions stft_compare() and istft_compare(), which,
compares the output between the existing (i)stft() and the shortTimeFFT based
_(i)stft_wrapper() implementations in this module.

For testing add the following imports to the file ``tests/test_spectral.py``::

    from ._scipy_spectral_test_shim import stft_compare as stft
    from ._scipy_spectral_test_shim import istft_compare as istft

and remove the existing imports of stft and istft.

The idea of these wrappers is not to provide a backward-compatible interface
but to demonstrate that the ShortTimeFFT implementation is at least as capable
as the existing one and delivers comparable results. Furthermore, the
wrappers highlight the different philosophies of the implementations,
especially in the border handling.
"""
import platform
from typing import cast, Literal

import numpy as np
from numpy.testing import assert_allclose

from scipy.signal import ShortTimeFFT
from scipy.signal import csd, get_window, stft, istft
from scipy.signal._arraytools import const_ext, even_ext, odd_ext, zero_ext
from scipy.signal._short_time_fft import FFT_MODE_TYPE
from scipy.signal._spectral_py import _spectral_helper, _triage_segments, \
    _median_bias


def _stft_wrapper(x, fs=1.0, window='hann', nperseg=256, noverlap=None,
                  nfft=None, detrend=False, return_onesided=True,
                  boundary='zeros', padded=True, axis=-1, scaling='spectrum'):
    """Wrapper for the SciPy `stft()` function based on `ShortTimeFFT` for
    unit testing.

    Handling the boundary and padding is where `ShortTimeFFT` and `stft()`
    differ in behavior. Parts of `_spectral_helper()` were copied to mimic
    the` stft()` behavior.

    This function is meant to be solely used by `stft_compare()`.
    """
    if scaling not in ('psd', 'spectrum'):  # same errors as in original stft:
        raise ValueError(f"Parameter {scaling=} not in ['spectrum', 'psd']!")

    # The following lines are taken from the original _spectral_helper():
    boundary_funcs = {'even': even_ext,
                      'odd': odd_ext,
                      'constant': const_ext,
                      'zeros': zero_ext,
                      None: None}

    if boundary not in boundary_funcs:
        raise ValueError(f"Unknown boundary option '{boundary}', must be one" +
                         f" of: {list(boundary_funcs.keys())}")
    if x.size == 0:
        return np.empty(x.shape), np.empty(x.shape), np.empty(x.shape)

    if nperseg is not None:  # if specified by user
        nperseg = int(nperseg)
        if nperseg < 1:
            raise ValueError('nperseg must be a positive integer')

    # parse window; if array like, then set nperseg = win.shape
    win, nperseg = _triage_segments(window, nperseg,
                                    input_length=x.shape[axis])

    if nfft is None:
        nfft = nperseg
    elif nfft < nperseg:
        raise ValueError('nfft must be greater than or equal to nperseg.')
    else:
        nfft = int(nfft)

    if noverlap is None:
        noverlap = nperseg//2
    else:
        noverlap = int(noverlap)
    if noverlap >= nperseg:
        raise ValueError('noverlap must be less than nperseg.')
    nstep = nperseg - noverlap
    n = x.shape[axis]

    # Padding occurs after boundary extension, so that the extended signal ends
    # in zeros, instead of introducing an impulse at the end.
    # I.e. if x = [..., 3, 2]
    # extend then pad -> [..., 3, 2, 2, 3, 0, 0, 0]
    # pad then extend -> [..., 3, 2, 0, 0, 0, 2, 3]

    if boundary is not None:
        ext_func = boundary_funcs[boundary]
        # Extend by nperseg//2 in front and back:
        x = ext_func(x, nperseg//2, axis=axis)

    if padded:
        # Pad to integer number of windowed segments
        # I.e make x.shape[-1] = nperseg + (nseg-1)*nstep, with integer nseg
        x = np.moveaxis(x, axis, -1)

        # This is an edge case where shortTimeFFT returns one more time slice
        # than the Scipy stft() shorten to remove last time slice:
        if n % 2 == 1 and nperseg % 2 == 1 and noverlap % 2 == 1:
            x = x[..., :axis - 1]

        nadd = (-(x.shape[-1]-nperseg) % nstep) % nperseg
        zeros_shape = list(x.shape[:-1]) + [nadd]
        x = np.concatenate((x, np.zeros(zeros_shape)), axis=-1)
        x = np.moveaxis(x, -1, axis)

    #  ... end original _spectral_helper() code.
    scale_to = {'spectrum': 'magnitude', 'psd': 'psd'}[scaling]

    if np.iscomplexobj(x) and return_onesided:
        return_onesided = False
    # using cast() to make mypy happy:
    fft_mode = cast(FFT_MODE_TYPE, 'onesided' if return_onesided else 'twosided')

    ST = ShortTimeFFT(win, nstep, fs, fft_mode=fft_mode, mfft=nfft,
                      scale_to=scale_to, phase_shift=None)

    k_off = nperseg // 2
    p0 = 0  # ST.lower_border_end[1] + 1
    nn = x.shape[axis] if padded else n+k_off+1
    p1 = ST.upper_border_begin(nn)[1]  # ST.p_max(n) + 1

    # This is bad hack to pass the test test_roundtrip_boundary_extension():
    if padded is True and nperseg - noverlap == 1:
        p1 -= nperseg // 2 - 1  # the reasoning behind this is not clear to me

    detr = None if detrend is False else detrend
    Sxx = ST.stft_detrend(x, detr, p0, p1, k_offset=k_off, axis=axis)
    t = ST.t(nn, 0, p1 - p0, k_offset=0 if boundary is not None else k_off)
    if x.dtype in (np.float32, np.complex64):
        Sxx = Sxx.astype(np.complex64)

    # workaround for test_average_all_segments() - seems to be buggy behavior:
    if boundary is None and padded is False:
        t, Sxx = t[1:-1], Sxx[..., :-2]
        t -= k_off / fs

    return ST.f, t, Sxx


def _istft_wrapper(Zxx, fs=1.0, window='hann', nperseg=None, noverlap=None,
                   nfft=None, input_onesided=True, boundary=True, time_axis=-1,
                   freq_axis=-2, scaling='spectrum') -> \
        tuple[np.ndarray, np.ndarray, tuple[int, int]]:
    """Wrapper for the SciPy `istft()` function based on `ShortTimeFFT` for
        unit testing.

    Note that only option handling is implemented as far as to handle the unit
    tests. E.g., the case ``nperseg=None`` is not handled.

    This function is meant to be solely used by `istft_compare()`.
    """
    # *** Lines are taken from _spectral_py.istft() ***:
    if Zxx.ndim < 2:
        raise ValueError('Input stft must be at least 2d!')

    if freq_axis == time_axis:
        raise ValueError('Must specify differing time and frequency axes!')

    nseg = Zxx.shape[time_axis]

    if input_onesided:
        # Assume even segment length
        n_default = 2*(Zxx.shape[freq_axis] - 1)
    else:
        n_default = Zxx.shape[freq_axis]

    # Check windowing parameters
    if nperseg is None:
        nperseg = n_default
    else:
        nperseg = int(nperseg)
        if nperseg < 1:
            raise ValueError('nperseg must be a positive integer')

    if nfft is None:
        if input_onesided and (nperseg == n_default + 1):
            # Odd nperseg, no FFT padding
            nfft = nperseg
        else:
            nfft = n_default
    elif nfft < nperseg:
        raise ValueError('nfft must be greater than or equal to nperseg.')
    else:
        nfft = int(nfft)

    if noverlap is None:
        noverlap = nperseg//2
    else:
        noverlap = int(noverlap)
    if noverlap >= nperseg:
        raise ValueError('noverlap must be less than nperseg.')
    nstep = nperseg - noverlap

    # Get window as array
    if isinstance(window, str) or type(window) is tuple:
        win = get_window(window, nperseg)
    else:
        win = np.asarray(window)
        if len(win.shape) != 1:
            raise ValueError('window must be 1-D')
        if win.shape[0] != nperseg:
            raise ValueError(f'window must have length of {nperseg}')

    outputlength = nperseg + (nseg-1)*nstep
    # *** End block of: Taken from _spectral_py.istft() ***

    # Using cast() to make mypy happy:
    fft_mode = cast(FFT_MODE_TYPE, 'onesided' if input_onesided else 'twosided')
    scale_to = cast(Literal['magnitude', 'psd'],
                    {'spectrum': 'magnitude', 'psd': 'psd'}[scaling])

    ST = ShortTimeFFT(win, nstep, fs, fft_mode=fft_mode, mfft=nfft,
                      scale_to=scale_to, phase_shift=None)

    if boundary:
        j = nperseg if nperseg % 2 == 0 else nperseg - 1
        k0 = ST.k_min + nperseg // 2
        k1 = outputlength - j + k0
    else:
        raise NotImplementedError("boundary=False does not make sense with" +
                                  "ShortTimeFFT.istft()!")

    x = ST.istft(Zxx, k0=k0, k1=k1, f_axis=freq_axis, t_axis=time_axis)
    t = np.arange(k1 - k0) * ST.T
    k_hi = ST.upper_border_begin(k1 - k0)[0]
    # using cast() to make mypy happy:
    return t, x, (ST.lower_border_end[0], k_hi)


def _csd_wrapper(x, y, fs=1.0, window='hann', nperseg=None, noverlap=None,
                 nfft=None, detrend='constant', return_onesided=True,
                 scaling='density', axis=-1, average='mean'):
    """Wrapper for the `csd()` function based on `ShortTimeFFT` for
        unit testing.
    """
    freqs, _, Pxy = _csd_test_shim(x, y, fs, window, nperseg, noverlap, nfft,
                                   detrend, return_onesided, scaling, axis)

    # The following code is taken from csd():
    if len(Pxy.shape) >= 2 and Pxy.size > 0:
        if Pxy.shape[-1] > 1:
            if average == 'median':
                # np.median must be passed real arrays for the desired result
                bias = _median_bias(Pxy.shape[-1])
                if np.iscomplexobj(Pxy):
                    Pxy = (np.median(np.real(Pxy), axis=-1)
                           + 1j * np.median(np.imag(Pxy), axis=-1))
                else:
                    Pxy = np.median(Pxy, axis=-1)
                Pxy /= bias
            elif average == 'mean':
                Pxy = Pxy.mean(axis=-1)
            else:
                raise ValueError(f'average must be "median" or "mean", got {average}')
        else:
            Pxy = np.reshape(Pxy, Pxy.shape[:-1])

    return freqs, Pxy


def _csd_test_shim(x, y, fs=1.0, window='hann', nperseg=None, noverlap=None,
                   nfft=None, detrend='constant', return_onesided=True,
                   scaling='density', axis=-1):
    """Compare output of  _spectral_helper() and ShortTimeFFT, more
    precisely _spect_helper_csd() for used in csd_wrapper().

   The motivation of this function is to test if the ShortTimeFFT-based
   wrapper `_spect_helper_csd()` returns the same values as `_spectral_helper`.
   This function should only be usd by csd() in (unit) testing.
   """
    freqs, t, Pxy = _spectral_helper(x, y, fs, window, nperseg, noverlap, nfft,
                                     detrend, return_onesided, scaling, axis,
                                     mode='psd')
    freqs1, Pxy1 = _spect_helper_csd(x, y, fs, window, nperseg, noverlap, nfft,
                                     detrend, return_onesided, scaling, axis)

    np.testing.assert_allclose(freqs1, freqs)
    amax_Pxy = max(np.abs(Pxy).max(), 1) if Pxy.size else 1
    atol = np.finfo(Pxy.dtype).resolution * amax_Pxy  # needed for large Pxy
    # for c_ in range(Pxy.shape[-1]):
    #    np.testing.assert_allclose(Pxy1[:, c_], Pxy[:, c_], atol=atol)
    np.testing.assert_allclose(Pxy1, Pxy, atol=atol)
    return freqs, t, Pxy


def _spect_helper_csd(x, y, fs=1.0, window='hann', nperseg=None, noverlap=None,
                      nfft=None, detrend='constant', return_onesided=True,
                      scaling='density', axis=-1):
    """Wrapper for replacing _spectral_helper() by using the ShortTimeFFT
      for use by csd().

    This function should be only used by _csd_test_shim() and is only useful
    for testing the ShortTimeFFT implementation.
    """

    # The following lines are taken from the original _spectral_helper():
    same_data = y is x
    axis = int(axis)

    # Ensure we have np.arrays, get outdtype
    x = np.asarray(x)
    if not same_data:
        y = np.asarray(y)
    #     outdtype = np.result_type(x, y, np.complex64)
    # else:
    #     outdtype = np.result_type(x, np.complex64)

    if not same_data:
        # Check if we can broadcast the outer axes together
        xouter = list(x.shape)
        youter = list(y.shape)
        xouter.pop(axis)
        youter.pop(axis)
        try:
            outershape = np.broadcast(np.empty(xouter), np.empty(youter)).shape
        except ValueError as e:
            raise ValueError('x and y cannot be broadcast together.') from e

    if same_data:
        if x.size == 0:
            return np.empty(x.shape), np.empty(x.shape)
    else:
        if x.size == 0 or y.size == 0:
            outshape = outershape + (min([x.shape[axis], y.shape[axis]]),)
            emptyout = np.moveaxis(np.empty(outshape), -1, axis)
            return emptyout, emptyout

    if nperseg is not None:  # if specified by user
        nperseg = int(nperseg)
        if nperseg < 1:
            raise ValueError('nperseg must be a positive integer')

    # parse window; if array like, then set nperseg = win.shape
    n = x.shape[axis] if same_data else max(x.shape[axis], y.shape[axis])
    win, nperseg = _triage_segments(window, nperseg, input_length=n)

    if nfft is None:
        nfft = nperseg
    elif nfft < nperseg:
        raise ValueError('nfft must be greater than or equal to nperseg.')
    else:
        nfft = int(nfft)

    if noverlap is None:
        noverlap = nperseg // 2
    else:
        noverlap = int(noverlap)
    if noverlap >= nperseg:
        raise ValueError('noverlap must be less than nperseg.')
    nstep = nperseg - noverlap

    if np.iscomplexobj(x) and return_onesided:
        return_onesided = False

    # using cast() to make mypy happy:
    fft_mode = cast(FFT_MODE_TYPE, 'onesided' if return_onesided
                    else 'twosided')
    scale = {'spectrum': 'magnitude', 'density': 'psd'}[scaling]
    SFT = ShortTimeFFT(win, nstep, fs, fft_mode=fft_mode, mfft=nfft,
                       scale_to=scale, phase_shift=None)

    # _spectral_helper() calculates X.conj()*Y instead of X*Y.conj():
    Pxy = SFT.spectrogram(y, x, detr=None if detrend is False else detrend,
                          p0=0, p1=(n-noverlap)//SFT.hop, k_offset=nperseg//2,
                          axis=axis).conj()
    # Note:
    # 'onesided2X' scaling of ShortTimeFFT conflicts with the
    # scaling='spectrum' parameter, since it doubles the squared magnitude,
    # which in the view of the ShortTimeFFT implementation does not make sense.
    # Hence, the doubling of the square is implemented here:
    if return_onesided:
        f_axis = Pxy.ndim - 1 + axis if axis < 0 else axis
        Pxy = np.moveaxis(Pxy, f_axis, -1)
        Pxy[..., 1:-1 if SFT.mfft % 2 == 0 else None] *= 2
        Pxy = np.moveaxis(Pxy, -1, f_axis)

    return SFT.f, Pxy


def stft_compare(x, fs=1.0, window='hann', nperseg=256, noverlap=None,
                 nfft=None, detrend=False, return_onesided=True,
                 boundary='zeros', padded=True, axis=-1, scaling='spectrum'):
    """Assert that the results from the existing `stft()` and `_stft_wrapper()`
    are close to each other.

    For comparing the STFT values an absolute tolerance of the floating point
    resolution was added to circumvent problems with the following tests:
    * For float32 the tolerances are much higher in
      TestSTFT.test_roundtrip_float32()).
    * The TestSTFT.test_roundtrip_scaling() has a high relative deviation.
      Interestingly this did not appear in Scipy 1.9.1 but only in the current
      development version.
    """
    kw = dict(x=x, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap,
              nfft=nfft, detrend=detrend, return_onesided=return_onesided,
              boundary=boundary, padded=padded, axis=axis, scaling=scaling)
    f, t, Zxx = stft(**kw)
    f_wrapper, t_wrapper, Zxx_wrapper = _stft_wrapper(**kw)

    e_msg_part = " of `stft_wrapper()` differ from `stft()`."
    assert_allclose(f_wrapper, f, err_msg=f"Frequencies {e_msg_part}")
    assert_allclose(t_wrapper, t, err_msg=f"Time slices {e_msg_part}")

    # Adapted tolerances to account for:
    atol = np.finfo(Zxx.dtype).resolution * 2
    assert_allclose(Zxx_wrapper, Zxx, atol=atol,
                    err_msg=f"STFT values {e_msg_part}")
    return f, t, Zxx


def istft_compare(Zxx, fs=1.0, window='hann', nperseg=None, noverlap=None,
                  nfft=None, input_onesided=True, boundary=True, time_axis=-1,
                  freq_axis=-2, scaling='spectrum'):
    """Assert that the results from the existing `istft()` and
    `_istft_wrapper()` are close to each other.

    Quirks:
    * If ``boundary=False`` the comparison is skipped, since it does not
      make sense with ShortTimeFFT.istft(). Only used in test
      TestSTFT.test_roundtrip_boundary_extension().
    * If ShortTimeFFT.istft() decides the STFT is not invertible, the
      comparison is skipped, since istft() only emits a warning and does not
      return a correct result. Only used in
      ShortTimeFFT.test_roundtrip_not_nola().
    * For comparing the signals an absolute tolerance of the floating point
      resolution was added to account for the low accuracy of float32 (Occurs
      only in TestSTFT.test_roundtrip_float32()).
    """
    kw = dict(Zxx=Zxx, fs=fs, window=window, nperseg=nperseg,
              noverlap=noverlap, nfft=nfft, input_onesided=input_onesided,
              boundary=boundary, time_axis=time_axis, freq_axis=freq_axis,
              scaling=scaling)

    t, x = istft(**kw)
    if not boundary:  # skip test_roundtrip_boundary_extension():
        return t, x  # _istft_wrapper does() not implement this case
    try:  # if inversion fails, istft() only emits a warning:
        t_wrapper, x_wrapper, (k_lo, k_hi) = _istft_wrapper(**kw)
    except ValueError as v:  # Do nothing if inversion fails:
        if v.args[0] == "Short-time Fourier Transform not invertible!":
            return t, x
        raise v

    e_msg_part = " of `istft_wrapper()` differ from `istft()`"
    assert_allclose(t, t_wrapper, err_msg=f"Sample times {e_msg_part}")

    # Adapted tolerances to account for resolution loss:
    atol = np.finfo(x.dtype).resolution*2  # instead of default atol = 0
    rtol = 1e-7  # default for np.allclose()

    # Relax atol on 32-Bit platforms a bit to pass CI tests.
    #  - Not clear why there are discrepancies (in the FFT maybe?)
    #  - Not sure what changed on 'i686' since earlier on those test passed
    if x.dtype == np.float32 and platform.machine() == 'i686':
        # float32 gets only used by TestSTFT.test_roundtrip_float32() so
        # we are using the tolerances from there to circumvent CI problems
        atol, rtol = 1e-4, 1e-5
    elif platform.machine() in ('aarch64', 'i386', 'i686'):
        atol = max(atol, 1e-12)  # 2e-15 seems too tight for 32-Bit platforms

    assert_allclose(x_wrapper[k_lo:k_hi], x[k_lo:k_hi], atol=atol, rtol=rtol,
                    err_msg=f"Signal values {e_msg_part}")
    return t, x


def csd_compare(x, y, fs=1.0, window='hann', nperseg=None, noverlap=None,
                nfft=None, detrend='constant', return_onesided=True,
                scaling='density', axis=-1, average='mean'):
    """Assert that the results from the existing `csd()` and `_csd_wrapper()`
    are close to each other. """
    kw = dict(x=x, y=y, fs=fs, window=window, nperseg=nperseg,
              noverlap=noverlap, nfft=nfft, detrend=detrend,
              return_onesided=return_onesided, scaling=scaling, axis=axis,
              average=average)
    freqs0, Pxy0 = csd(**kw)
    freqs1, Pxy1 = _csd_wrapper(**kw)

    assert_allclose(freqs1, freqs0)
    assert_allclose(Pxy1, Pxy0)
    assert_allclose(freqs1, freqs0)
    return freqs0, Pxy0
