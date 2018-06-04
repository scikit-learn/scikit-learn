# Copyright (c) 2006-2012 Filip Wasilewski <http://en.ig.ma/>
# Copyright (c) 2012-2016 The PyWavelets Developers
#                         <https://github.com/PyWavelets/pywt>
# See COPYING for license details.

"""
2D and nD Discrete Wavelet Transforms and Inverse Discrete Wavelet Transforms.
"""

from __future__ import division, print_function, absolute_import

__all__ = ['dwt2', 'idwt2', 'dwtn', 'idwtn']

from itertools import product

import numpy as np

from ._extensions._pywt import Wavelet, Modes
from ._extensions._dwt import dwt_axis, idwt_axis


def dwt2(data, wavelet, mode='symmetric', axes=(-2, -1)):
    """
    2D Discrete Wavelet Transform.

    Parameters
    ----------
    data : array_like
        2D array with input data
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    axes : 2-tuple of ints, optional
        Axes over which to compute the DWT. Repeated elements mean the DWT will
        be performed multiple times along these axes.

    Returns
    -------
    (cA, (cH, cV, cD)) : tuple
        Approximation, horizontal detail, vertical detail and diagonal
        detail coefficients respectively.  Horizontal refers to array axis 0
        (or ``axes[0]`` for user-specified ``axes``).

    Examples
    --------
    >>> import numpy as np
    >>> import pywt
    >>> data = np.ones((4,4), dtype=np.float64)
    >>> coeffs = pywt.dwt2(data, 'haar')
    >>> cA, (cH, cV, cD) = coeffs
    >>> cA
    array([[ 2.,  2.],
           [ 2.,  2.]])
    >>> cV
    array([[ 0.,  0.],
           [ 0.,  0.]])

    """
    axes = tuple(axes)
    data = np.asarray(data)
    if len(axes) != 2:
        raise ValueError("Expected 2 axes")
    if data.ndim < len(np.unique(axes)):
        raise ValueError("Input array has fewer dimensions than the specified "
                         "axes")

    coefs = dwtn(data, wavelet, mode, axes)
    return coefs['aa'], (coefs['da'], coefs['ad'], coefs['dd'])


def idwt2(coeffs, wavelet, mode='symmetric', axes=(-2, -1)):
    """
    2-D Inverse Discrete Wavelet Transform.

    Reconstructs data from coefficient arrays.

    Parameters
    ----------
    coeffs : tuple
        (cA, (cH, cV, cD)) A tuple with approximation coefficients and three
        details coefficients 2D arrays like from `dwt2`. If any of these
        components are set to ``None``, it will be treated as zeros.
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    axes : 2-tuple of ints, optional
        Axes over which to compute the IDWT. Repeated elements mean the IDWT
        will be performed multiple times along these axes.

    Examples
    --------
    >>> import numpy as np
    >>> import pywt
    >>> data = np.array([[1,2], [3,4]], dtype=np.float64)
    >>> coeffs = pywt.dwt2(data, 'haar')
    >>> pywt.idwt2(coeffs, 'haar')
    array([[ 1.,  2.],
           [ 3.,  4.]])

    """
    # L -low-pass data, H - high-pass data
    LL, (HL, LH, HH) = coeffs
    axes = tuple(axes)
    if len(axes) != 2:
        raise ValueError("Expected 2 axes")

    coeffs = {'aa': LL, 'da': HL, 'ad': LH, 'dd': HH}
    return idwtn(coeffs, wavelet, mode, axes)


def dwtn(data, wavelet, mode='symmetric', axes=None):
    """
    Single-level n-dimensional Discrete Wavelet Transform.

    Parameters
    ----------
    data : array_like
        n-dimensional array with input data.
    wavelet : Wavelet object or name string
        Wavelet to use.
    mode : str, optional
        Signal extension mode, see `Modes`.  Default is 'symmetric'.
    axes : sequence of ints, optional
        Axes over which to compute the DWT. Repeated elements mean the DWT will
        be performed multiple times along these axes. A value of ``None`` (the
        default) selects all axes.

        Axes may be repeated, but information about the original size may be
        lost if it is not divisible by ``2 ** nrepeats``. The reconstruction
        will be larger, with additional values derived according to the
        ``mode`` parameter. ``pywt.wavedecn`` should be used for multilevel
        decomposition.

    Returns
    -------
    coeffs : dict
        Results are arranged in a dictionary, where key specifies
        the transform type on each dimension and value is a n-dimensional
        coefficients array.

        For example, for a 2D case the result will look something like this::

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
        real = dwtn(data.real, wavelet, mode, axes)
        imag = dwtn(data.imag, wavelet, mode, axes)
        return dict((k, real[k] + 1j * imag[k]) for k in real.keys())

    if data.dtype == np.dtype('object'):
        raise TypeError("Input must be a numeric array-like")
    if data.ndim < 1:
        raise ValueError("Input data must be at least 1D")

    if axes is None:
        axes = range(data.ndim)
    axes = (a + data.ndim if a < 0 else a for a in axes)

    mode = Modes.from_object(mode)
    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    coeffs = [('', data)]
    for axis in axes:
        new_coeffs = []
        for subband, x in coeffs:
            cA, cD = dwt_axis(x, wavelet, mode, axis)
            new_coeffs.extend([(subband + 'a', cA),
                               (subband + 'd', cD)])
        coeffs = new_coeffs
    return dict(coeffs)


def _fix_coeffs(coeffs):
    missing_keys = [k for k, v in coeffs.items() if
                    v is None]
    if missing_keys:
        raise ValueError(
            "The following detail coefficients were set to None: "
            "{}.".format(missing_keys))

    invalid_keys = [k for k, v in coeffs.items() if
                    not set(k) <= set('ad')]
    if invalid_keys:
        raise ValueError(
            "The following invalid keys were found in the detail "
            "coefficient dictionary: {}.".format(invalid_keys))

    key_lengths = [len(k) for k in coeffs.keys()]
    if len(np.unique(key_lengths)) > 1:
        raise ValueError(
            "All detail coefficient names must have equal length.")

    return dict((k, np.asarray(v)) for k, v in coeffs.items())


def idwtn(coeffs, wavelet, mode='symmetric', axes=None):
    """
    Single-level n-dimensional Inverse Discrete Wavelet Transform.

    Parameters
    ----------
    coeffs: dict
        Dictionary as in output of `dwtn`. Missing or ``None`` items
        will be treated as zeros.
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode used in the decomposition,
        see Modes (default: 'symmetric').
    axes : sequence of ints, optional
        Axes over which to compute the IDWT. Repeated elements mean the IDWT
        will be performed multiple times along these axes. A value of ``None``
        (the default) selects all axes.

        For the most accurate reconstruction, the axes should be provided in
        the same order as they were provided to ``dwtn``.

    Returns
    -------
    data: ndarray
        Original signal reconstructed from input data.

    """
    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)
    mode = Modes.from_object(mode)

    # drop the keys corresponding to value = None
    coeffs = dict((k, v) for k, v in coeffs.items() if v is not None)

    # Raise error for invalid key combinations
    coeffs = _fix_coeffs(coeffs)

    if any(np.iscomplexobj(v) for v in coeffs.values()):
        real_coeffs = dict((k, v.real) for k, v in coeffs.items())
        imag_coeffs = dict((k, v.imag) for k, v in coeffs.items())
        return (idwtn(real_coeffs, wavelet, mode, axes) +
                1j * idwtn(imag_coeffs, wavelet, mode, axes))

    # key length matches the number of axes transformed
    ndim_transform = max(len(key) for key in coeffs.keys())

    try:
        coeff_shapes = (v.shape for k, v in coeffs.items()
                        if v is not None and len(k) == ndim_transform)
        coeff_shape = next(coeff_shapes)
    except StopIteration:
        raise ValueError("`coeffs` must contain at least one non-null wavelet "
                         "band")
    if any(s != coeff_shape for s in coeff_shapes):
        raise ValueError("`coeffs` must all be of equal size (or None)")

    if axes is None:
        axes = range(ndim_transform)
        ndim = ndim_transform
    else:
        ndim = len(coeff_shape)
    axes = (a + ndim if a < 0 else a for a in axes)

    for key_length, axis in reversed(list(enumerate(axes))):
        if axis < 0 or axis >= ndim:
            raise ValueError("Axis greater than data dimensions")

        new_coeffs = {}
        new_keys = [''.join(coef) for coef in product('ad', repeat=key_length)]

        for key in new_keys:
            L = coeffs.get(key + 'a', None)
            H = coeffs.get(key + 'd', None)

            new_coeffs[key] = idwt_axis(L, H, wavelet, mode, axis)
        coeffs = new_coeffs

    return coeffs['']
