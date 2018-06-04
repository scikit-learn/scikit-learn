from ._extensions._swt import swt_max_level, swt as _swt, swt_axis as _swt_axis
from ._extensions._pywt import Wavelet, _check_dtype

import warnings
import numpy as np
from ._extensions._swt import swt_axis as _swt_axis
from ._dwt import idwt
from ._multidim import idwt2

__all__ = ["swt", "swt_max_level", 'iswt', 'swt2', 'iswt2', 'swtn']


def swt(data, wavelet, level=None, start_level=0, axis=-1):
    """
    Multilevel 1D stationary wavelet transform.

    Parameters
    ----------
    data :
        Input signal
    wavelet :
        Wavelet to use (Wavelet object or name)
    level : int, optional
        The number of decomposition steps to perform.
    start_level : int, optional
        The level at which the decomposition will begin (it allows one to
        skip a given number of transform steps and compute
        coefficients starting from start_level) (default: 0)
    axis: int, optional
        Axis over which to compute the SWT. If not given, the
        last axis is used.

    Returns
    -------
    coeffs : list
        List of approximation and details coefficients pairs in order
        similar to wavedec function::

            [(cAn, cDn), ..., (cA2, cD2), (cA1, cD1)]

        where n equals input parameter ``level``.

        If ``start_level = m`` is given, then the beginning m steps are
        skipped::

            [(cAm+n, cDm+n), ..., (cAm+1, cDm+1), (cAm, cDm)]

    """
    if np.iscomplexobj(data):
        data = np.asarray(data)
        coeffs_real = swt(data.real, wavelet, level, start_level)
        coeffs_imag = swt(data.imag, wavelet, level, start_level)
        coeffs_cplx = []
        for (cA_r, cD_r), (cA_i, cD_i) in zip(coeffs_real, coeffs_imag):
            coeffs_cplx.append((cA_r + 1j*cA_i, cD_r + 1j*cD_i))
        return coeffs_cplx

    # accept array_like input; make a copy to ensure a contiguous array
    dt = _check_dtype(data)
    data = np.array(data, dtype=dt)
    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)
    if level is None:
        level = swt_max_level(len(data))

    if axis < 0:
        axis = axis + data.ndim
    if not 0 <= axis < data.ndim:
        raise ValueError("Axis greater than data dimensions")

    if data.ndim == 1:
        ret = _swt(data, wavelet, level, start_level)
    else:
        ret = _swt_axis(data, wavelet, level, start_level, axis)
    return [(np.asarray(cA), np.asarray(cD)) for cA, cD in ret]


def iswt(coeffs, wavelet):
    """
    Multilevel 1D inverse discrete stationary wavelet transform.

    Parameters
    ----------
    coeffs : array_like
        Coefficients list of tuples::

            [(cAn, cDn), ..., (cA2, cD2), (cA1, cD1)]

        where cA is approximation, cD is details.  Index 1 corresponds to
        ``start_level`` from ``pywt.swt``.
    wavelet : Wavelet object or name string
        Wavelet to use

    Returns
    -------
    1D array of reconstructed data.

    Examples
    --------
    >>> import pywt
    >>> coeffs = pywt.swt([1,2,3,4,5,6,7,8], 'db2', level=2)
    >>> pywt.iswt(coeffs, 'db2')
    array([ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.])
    """

    output = coeffs[0][0].copy()  # Avoid modification of input data

    # num_levels, equivalent to the decomposition level, n
    num_levels = len(coeffs)
    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)
    for j in range(num_levels, 0, -1):
        step_size = int(pow(2, j-1))
        last_index = step_size
        _, cD = coeffs[num_levels - j]
        for first in range(last_index):  # 0 to last_index - 1

            # Getting the indices that we will transform
            indices = np.arange(first, len(cD), step_size)

            # select the even indices
            even_indices = indices[0::2]
            # select the odd indices
            odd_indices = indices[1::2]

            # perform the inverse dwt on the selected indices,
            # making sure to use periodic boundary conditions
            x1 = idwt(output[even_indices], cD[even_indices],
                      wavelet, 'periodization')
            x2 = idwt(output[odd_indices], cD[odd_indices],
                      wavelet, 'periodization')

            # perform a circular shift right
            x2 = np.roll(x2, 1)

            # average and insert into the correct indices
            output[indices] = (x1 + x2)/2.

    return output


def swt2(data, wavelet, level, start_level=0, axes=(-2, -1)):
    """
    Multilevel 2D stationary wavelet transform.

    Parameters
    ----------
    data : array_like
        2D array with input data
    wavelet : Wavelet object or name string
        Wavelet to use
    level : int
        The number of decomposition steps to perform.
    start_level : int, optional
        The level at which the decomposition will start (default: 0)
    axes : 2-tuple of ints, optional
        Axes over which to compute the SWT. Repeated elements are not allowed.

    Returns
    -------
    coeffs : list
        Approximation and details coefficients::

            [
                (cA_m,
                    (cH_m, cV_m, cD_m)
                ),
                (cA_m+1,
                    (cH_m+1, cV_m+1, cD_m+1)
                ),
                ...,
                (cA_m+level,
                    (cH_m+level, cV_m+level, cD_m+level)
                )
            ]

        where cA is approximation, cH is horizontal details, cV is
        vertical details, cD is diagonal details and m is ``start_level``.

    """
    axes = tuple(axes)
    data = np.asarray(data)
    if len(axes) != 2:
        raise ValueError("Expected 2 axes")
    if len(axes) != len(set(axes)):
        raise ValueError("The axes passed to swt2 must be unique.")
    if data.ndim < len(np.unique(axes)):
        raise ValueError("Input array has fewer dimensions than the specified "
                         "axes")

    coefs = swtn(data, wavelet, level, start_level, axes)
    ret = []
    for c in coefs:
        ret.append((c['aa'], (c['da'], c['ad'], c['dd'])))

    warnings.warn(
        "For consistency with the rest of PyWavelets, the order of the list "
        "returned by swt2 will be reversed in the next release. "
        "In other words, the levels will be sorted in descending rather than "
        "ascending order.", FutureWarning)
    # reverse order for backwards compatiblity
    ret.reverse()
    return ret


def iswt2(coeffs, wavelet):
    """
    Multilevel 2D inverse discrete stationary wavelet transform.

    Parameters
    ----------
    coeffs : list
        Approximation and details coefficients::

            [
                (cA_1,
                    (cH_1, cV_1, cD_1)
                ),
                (cA_2,
                    (cH_2, cV_2, cD_2)
                ),
                ...,
                (cA_n
                    (cH_n, cV_n, cD_n)
                )
            ]

        where cA is approximation, cH is horizontal details, cV is
        vertical details, cD is diagonal details and n is the number of
        levels.  Index 1 corresponds to ``start_level`` from ``pywt.swt2``.
    wavelet : Wavelet object or name string
        Wavelet to use

    Returns
    -------
    2D array of reconstructed data.

    Examples
    --------
    >>> import pywt
    >>> coeffs = pywt.swt2([[1,2,3,4],[5,6,7,8],
    ...                     [9,10,11,12],[13,14,15,16]],
    ...                    'db1', level=2)
    >>> pywt.iswt2(coeffs, 'db1')
    array([[  1.,   2.,   3.,   4.],
           [  5.,   6.,   7.,   8.],
           [  9.,  10.,  11.,  12.],
           [ 13.,  14.,  15.,  16.]])

    """

    output = coeffs[-1][0].copy()  # Avoid modification of input data

    # num_levels, equivalent to the decomposition level, n
    num_levels = len(coeffs)
    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)
    for j in range(num_levels, 0, -1):
        step_size = int(pow(2, j-1))
        last_index = step_size
        _, (cH, cV, cD) = coeffs[j-1]
        # We are going to assume cH, cV, and cD are square and of equal size
        if (cH.shape != cV.shape) or (cH.shape != cD.shape) or (
                cH.shape[0] != cH.shape[1]):
            raise RuntimeError(
                "Mismatch in shape of intermediate coefficient arrays")
        for first_h in range(last_index):  # 0 to last_index - 1
            for first_w in range(last_index):  # 0 to last_index - 1
                # Getting the indices that we will transform
                indices_h = slice(first_h, cH.shape[0], step_size)
                indices_w = slice(first_w, cH.shape[1], step_size)

                even_idx_h = slice(first_h, cH.shape[0], 2*step_size)
                even_idx_w = slice(first_w, cH.shape[1], 2*step_size)
                odd_idx_h = slice(first_h + step_size, cH.shape[0], 2*step_size)
                odd_idx_w = slice(first_w + step_size, cH.shape[1], 2*step_size)

                # perform the inverse dwt on the selected indices,
                # making sure to use periodic boundary conditions
                x1 = idwt2((output[even_idx_h, even_idx_w],
                           (cH[even_idx_h, even_idx_w],
                            cV[even_idx_h, even_idx_w],
                            cD[even_idx_h, even_idx_w])),
                           wavelet, 'periodization')
                x2 = idwt2((output[even_idx_h, odd_idx_w],
                           (cH[even_idx_h, odd_idx_w],
                            cV[even_idx_h, odd_idx_w],
                            cD[even_idx_h, odd_idx_w])),
                           wavelet, 'periodization')
                x3 = idwt2((output[odd_idx_h, even_idx_w],
                           (cH[odd_idx_h, even_idx_w],
                            cV[odd_idx_h, even_idx_w],
                            cD[odd_idx_h, even_idx_w])),
                           wavelet, 'periodization')
                x4 = idwt2((output[odd_idx_h, odd_idx_w],
                           (cH[odd_idx_h, odd_idx_w],
                            cV[odd_idx_h, odd_idx_w],
                            cD[odd_idx_h, odd_idx_w])),
                           wavelet, 'periodization')

                # perform a circular shifts
                x2 = np.roll(x2, 1, axis=1)
                x3 = np.roll(x3, 1, axis=0)
                x4 = np.roll(x4, 1, axis=0)
                x4 = np.roll(x4, 1, axis=1)
                output[indices_h, indices_w] = (x1 + x2 + x3 + x4) / 4

    warnings.warn(
        "For consistency with the rest of PyWavelets, the order of levels in "
        "coeffs expected by iswt2 will be reversed in the next release. "
        "In other words, the levels will be sorted in descending rather than "
        "ascending order.", FutureWarning)

    return output


def swtn(data, wavelet, level, start_level=0, axes=None):
    """
    n-dimensional stationary wavelet transform.

    Parameters
    ----------
    data : array_like
        n-dimensional array with input data.
    wavelet : Wavelet object or name string
        Wavelet to use.
    level : int
        The number of decomposition steps to perform.
    start_level : int, optional
        The level at which the decomposition will start (default: 0)
    axes : sequence of ints, optional
        Axes over which to compute the SWT. A value of ``None`` (the
        default) selects all axes. Axes may not be repeated.

    Returns
    -------
    [{coeffs_level_n}, ..., {coeffs_level_1}]: list of dict
        Results for each level are arranged in a dictionary, where the key
        specifies the transform type on each dimension and value is a
        n-dimensional coefficients array.

        For example, for a 2D case the result at a given level will look
        something like this::

            {'aa': <coeffs>  # A(LL) - approx. on 1st dim, approx. on 2nd dim
             'ad': <coeffs>  # V(LH) - approx. on 1st dim, det. on 2nd dim
             'da': <coeffs>  # H(HL) - det. on 1st dim, approx. on 2nd dim
             'dd': <coeffs>  # D(HH) - det. on 1st dim, det. on 2nd dim
            }

        For user-specified ``axes``, the order of the characters in the
        dictionary keys map to the specified ``axes``.

    """
    data = np.asarray(data)
    if np.iscomplexobj(data):
        real = swtn(data.real, wavelet, level, start_level, axes)
        imag = swtn(data.imag, wavelet, level, start_level, axes)
        cplx = []
        for rdict, idict in zip(real, imag):
            cplx.append(
                dict((k, rdict[k] + 1j * idict[k]) for k in rdict.keys()))
        return cplx

    if data.dtype == np.dtype('object'):
        raise TypeError("Input must be a numeric array-like")
    if data.ndim < 1:
        raise ValueError("Input data must be at least 1D")

    if axes is None:
        axes = range(data.ndim)
    axes = [a + data.ndim if a < 0 else a for a in axes]
    if len(axes) != len(set(axes)):
        raise ValueError("The axes passed to swtn must be unique.")
    num_axes = len(axes)

    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    ret = []
    for i in range(start_level, start_level + level):
        coeffs = [('', data)]
        for axis in axes:
            new_coeffs = []
            for subband, x in coeffs:
                cA, cD = _swt_axis(x, wavelet, level=1, start_level=i,
                                   axis=axis)[0]
                new_coeffs.extend([(subband + 'a', cA),
                                   (subband + 'd', cD)])
            coeffs = new_coeffs

        coeffs = dict(coeffs)
        ret.append(coeffs)

        # data for the next level is the approximation coeffs from this level
        data = coeffs['a' * num_axes]

    ret.reverse()
    return ret
