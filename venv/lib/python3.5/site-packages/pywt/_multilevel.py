# Copyright (c) 2006-2012 Filip Wasilewski <http://en.ig.ma/>
# Copyright (c) 2012-2016 The PyWavelets Developers
#                         <https://github.com/PyWavelets/pywt>
# See COPYING for license details.

"""
Multilevel 1D and 2D Discrete Wavelet Transform
and Inverse Discrete Wavelet Transform.
"""

from __future__ import division, print_function, absolute_import

import warnings
from copy import copy
import numpy as np

from ._extensions._pywt import Wavelet
from ._extensions._dwt import dwt_max_level
from ._dwt import dwt, idwt
from ._multidim import dwt2, idwt2, dwtn, idwtn, _fix_coeffs

__all__ = ['wavedec', 'waverec', 'wavedec2', 'waverec2', 'wavedecn',
           'waverecn', 'coeffs_to_array', 'array_to_coeffs']


def _check_level(size, dec_len, level):
    """
    Set the default decomposition level or check if requested level is valid.
    """
    if level is None:
        level = dwt_max_level(size, dec_len)
    elif level < 0:
        raise ValueError(
            "Level value of %d is too low . Minimum level is 0." % level)
    else:
        max_level = dwt_max_level(size, dec_len)
        if level > max_level:
            raise ValueError(
                "Level value of %d is too high.  Maximum allowed is %d." % (
                    level, max_level))
    return level


def wavedec(data, wavelet, mode='symmetric', level=None, axis=-1):
    """
    Multilevel 1D Discrete Wavelet Transform of data.

    Parameters
    ----------
    data: array_like
        Input data
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    level : int, optional
        Decomposition level (must be >= 0). If level is None (default) then it
        will be calculated using the ``dwt_max_level`` function.
    axis: int, optional
        Axis over which to compute the DWT. If not given, the
        last axis is used.

    Returns
    -------
    [cA_n, cD_n, cD_n-1, ..., cD2, cD1] : list
        Ordered list of coefficients arrays
        where `n` denotes the level of decomposition. The first element
        (`cA_n`) of the result is approximation coefficients array and the
        following elements (`cD_n` - `cD_1`) are details coefficients arrays.

    Examples
    --------
    >>> from pywt import wavedec
    >>> coeffs = wavedec([1,2,3,4,5,6,7,8], 'db1', level=2)
    >>> cA2, cD2, cD1 = coeffs
    >>> cD1
    array([-0.70710678, -0.70710678, -0.70710678, -0.70710678])
    >>> cD2
    array([-2., -2.])
    >>> cA2
    array([  5.,  13.])

    """
    data = np.asarray(data)

    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    try:
        axes_shape = data.shape[axis]
    except IndexError:
        raise ValueError("Axis greater than data dimensions")
    level = _check_level(axes_shape, wavelet.dec_len, level)

    coeffs_list = []

    a = data
    for i in range(level):
        a, d = dwt(a, wavelet, mode, axis)
        coeffs_list.append(d)

    coeffs_list.append(a)
    coeffs_list.reverse()

    return coeffs_list


def waverec(coeffs, wavelet, mode='symmetric', axis=-1):
    """
    Multilevel 1D Inverse Discrete Wavelet Transform.

    Parameters
    ----------
    coeffs : array_like
        Coefficients list [cAn, cDn, cDn-1, ..., cD2, cD1]
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    axis: int, optional
        Axis over which to compute the inverse DWT. If not given, the
        last axis is used.

    Examples
    --------
    >>> import pywt
    >>> coeffs = pywt.wavedec([1,2,3,4,5,6,7,8], 'db1', level=2)
    >>> pywt.waverec(coeffs, 'db1')
    array([ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.])
    """

    if not isinstance(coeffs, (list, tuple)):
        raise ValueError("Expected sequence of coefficient arrays.")

    if len(coeffs) < 1:
        raise ValueError(
            "Coefficient list too short (minimum 1 arrays required).")
    elif len(coeffs) == 1:
        # level 0 transform (just returns the approximation coefficients)
        return coeffs[0]

    a, ds = coeffs[0], coeffs[1:]

    for d in ds:
        if (a is not None) and (d is not None):
            try:
                if a.shape[axis] == d.shape[axis] + 1:
                    a = a[[slice(s) for s in d.shape]]
                elif a.shape[axis] != d.shape[axis]:
                    raise ValueError("coefficient shape mismatch")
            except IndexError:
                raise ValueError("Axis greater than coefficient dimensions")
        a = idwt(a, d, wavelet, mode, axis)

    return a


def wavedec2(data, wavelet, mode='symmetric', level=None, axes=(-2, -1)):
    """
    Multilevel 2D Discrete Wavelet Transform.

    Parameters
    ----------
    data : ndarray
        2D input data
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    level : int, optional
        Decomposition level (must be >= 0). If level is None (default) then it
        will be calculated using the ``dwt_max_level`` function.
    axes : 2-tuple of ints, optional
        Axes over which to compute the DWT. Repeated elements are not allowed.

    Returns
    -------
    [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)] : list
        Coefficients list.  For user-specified ``axes``, ``cH*``
        corresponds to ``axes[0]`` while ``cV*`` corresponds to ``axes[1]``.
        The first element returned is the approximation coefficients for the
        nth level of decomposition.  Remaining elements are tuples of detail
        coefficients in descending order of decomposition level.
        (i.e. cH1 are the horizontal detail coefficients at the first level)

    Examples
    --------
    >>> import pywt
    >>> import numpy as np
    >>> coeffs = pywt.wavedec2(np.ones((4,4)), 'db1')
    >>> # Levels:
    >>> len(coeffs)-1
    2
    >>> pywt.waverec2(coeffs, 'db1')
    array([[ 1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.]])
    """
    data = np.asarray(data)
    if data.ndim < 2:
        raise ValueError("Expected input data to have at least 2 dimensions.")

    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    axes = tuple(axes)
    if len(axes) != 2:
        raise ValueError("Expected 2 axes")
    if len(axes) != len(set(axes)):
        raise ValueError("The axes passed to wavedec2 must be unique.")
    try:
        axes_sizes = [data.shape[ax] for ax in axes]
    except IndexError:
        raise ValueError("Axis greater than data dimensions")
    level = _check_level(min(axes_sizes), wavelet.dec_len, level)

    coeffs_list = []

    a = data
    for i in range(level):
        a, ds = dwt2(a, wavelet, mode, axes)
        coeffs_list.append(ds)

    coeffs_list.append(a)
    coeffs_list.reverse()

    return coeffs_list


def waverec2(coeffs, wavelet, mode='symmetric', axes=(-2, -1)):
    """
    Multilevel 2D Inverse Discrete Wavelet Transform.

    coeffs : list or tuple
        Coefficients list [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)]
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    axes : 2-tuple of ints, optional
        Axes over which to compute the IDWT. Repeated elements are not allowed.

    Returns
    -------
    2D array of reconstructed data.

    Examples
    --------
    >>> import pywt
    >>> import numpy as np
    >>> coeffs = pywt.wavedec2(np.ones((4,4)), 'db1')
    >>> # Levels:
    >>> len(coeffs)-1
    2
    >>> pywt.waverec2(coeffs, 'db1')
    array([[ 1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.]])
    """

    if not isinstance(coeffs, (list, tuple)):
        raise ValueError("Expected sequence of coefficient arrays.")

    if len(axes) != len(set(axes)):
        raise ValueError("The axes passed to waverec2 must be unique.")

    if len(coeffs) < 1:
        raise ValueError(
            "Coefficient list too short (minimum 1 array required).")
    elif len(coeffs) == 1:
        # level 0 transform (just returns the approximation coefficients)
        return coeffs[0]

    a, ds = coeffs[0], coeffs[1:]
    a = np.asarray(a)

    for d in ds:
        d = tuple(np.asarray(coeff) if coeff is not None else None
                  for coeff in d)
        d_shapes = (coeff.shape for coeff in d if coeff is not None)
        try:
            d_shape = next(d_shapes)
        except StopIteration:
            idxs = slice(None), slice(None)
        else:
            if not all(s == d_shape for s in d_shapes):
                raise ValueError("All detail shapes must be the same length.")
            idxs = tuple(slice(None, -1 if a_len == d_len + 1 else None)
                         for a_len, d_len in zip(a.shape, d_shape))
        a = idwt2((a[idxs], d), wavelet, mode, axes)

    return a


def wavedecn(data, wavelet, mode='symmetric', level=None, axes=None):
    """
    Multilevel nD Discrete Wavelet Transform.

    Parameters
    ----------
    data : ndarray
        nD input data
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    level : int, optional
        Decomposition level (must be >= 0). If level is None (default) then it
        will be calculated using the ``dwt_max_level`` function.
    axes : sequence of ints, optional
        Axes over which to compute the DWT. Axes may not be repeated. The
        default is ``None``, which means transform all axes
        (``axes = range(data.ndim)``).

    Returns
    -------
    [cAn, {details_level_n}, ... {details_level_1}] : list
        Coefficients list.  Coefficients are listed in descending order of
        decomposition level.  ``cAn`` are the approximation coefficients at
        level ``n``.  Each ``details_level_i`` element is a dictionary
        containing detail coefficients at level `i` of the decomposition.  As
        a concrete example, a 3D decomposition would have the following set of
        keys in each ``details_level_i`` dictionary::

            {'aad', 'ada', 'daa', 'add', 'dad', 'dda', 'ddd'}

        where the order of the characters in each key map to the specified
        ``axes``.

    Examples
    --------
    >>> import numpy as np
    >>> from pywt import wavedecn, waverecn
    >>> coeffs = wavedecn(np.ones((4, 4, 4)), 'db1')
    >>> # Levels:
    >>> len(coeffs)-1
    2
    >>> waverecn(coeffs, 'db1')  # doctest: +NORMALIZE_WHITESPACE
    array([[[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]],
           [[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]],
           [[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]],
           [[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]]])

    """
    data = np.asarray(data)

    if len(data.shape) < 1:
        raise ValueError("Expected at least 1D input data.")

    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    if np.isscalar(axes):
        axes = (axes, )
    if axes is None:
        axes = range(data.ndim)
    else:
        axes = tuple(axes)
    if len(axes) != len(set(axes)):
        raise ValueError("The axes passed to wavedecn must be unique.")
    ndim_transform = len(axes)
    try:
        axes_shapes = [data.shape[ax] for ax in axes]
    except IndexError:
        raise ValueError("Axis greater than data dimensions")
    level = _check_level(min(axes_shapes), wavelet.dec_len, level)

    coeffs_list = []

    a = data
    for i in range(level):
        coeffs = dwtn(a, wavelet, mode, axes)
        a = coeffs.pop('a' * ndim_transform)
        coeffs_list.append(coeffs)

    coeffs_list.append(a)
    coeffs_list.reverse()

    return coeffs_list


def _match_coeff_dims(a_coeff, d_coeff_dict):
    # For each axis, compare the approximation coeff shape to one of the
    # stored detail coeffs and truncate the last element along the axis
    # if necessary.
    if a_coeff is None:
        return None
    if not d_coeff_dict:
        return a_coeff
    d_coeff = d_coeff_dict[next(iter(d_coeff_dict))]
    size_diffs = np.subtract(a_coeff.shape, d_coeff.shape)
    if np.any((size_diffs < 0) | (size_diffs > 1)):
        raise ValueError("incompatible coefficient array sizes")
    return a_coeff[[slice(s) for s in d_coeff.shape]]


def waverecn(coeffs, wavelet, mode='symmetric', axes=None):
    """
    Multilevel nD Inverse Discrete Wavelet Transform.

    coeffs : array_like
        Coefficients list [cAn, {details_level_n}, ... {details_level_1}]
    wavelet : Wavelet object or name string
        Wavelet to use
    mode : str, optional
        Signal extension mode, see Modes (default: 'symmetric')
    axes : sequence of ints, optional
        Axes over which to compute the IDWT.  Axes may not be repeated.

    Returns
    -------
    nD array of reconstructed data.

    Examples
    --------
    >>> import numpy as np
    >>> from pywt import wavedecn, waverecn
    >>> coeffs = wavedecn(np.ones((4, 4, 4)), 'db1')
    >>> # Levels:
    >>> len(coeffs)-1
    2
    >>> waverecn(coeffs, 'db1')  # doctest: +NORMALIZE_WHITESPACE
    array([[[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]],
           [[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]],
           [[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]],
           [[ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.],
            [ 1.,  1.,  1.,  1.]]])

    """
    if len(coeffs) < 1:
        raise ValueError(
            "Coefficient list too short (minimum 1 array required).")

    a, ds = coeffs[0], coeffs[1:]

    # Raise error for invalid key combinations
    ds = list(map(_fix_coeffs, ds))

    if not ds:
        # level 0 transform (just returns the approximation coefficients)
        return coeffs[0]
    if a is None and not any(ds):
        raise ValueError("At least one coefficient must contain a valid value.")

    coeff_ndims = []
    if a is not None:
        a = np.asarray(a)
        coeff_ndims.append(a.ndim)
    for d in ds:
        coeff_ndims += [v.ndim for k, v in d.items()]

    # test that all coefficients have a matching number of dimensions
    unique_coeff_ndims = np.unique(coeff_ndims)
    if len(unique_coeff_ndims) == 1:
        ndim = unique_coeff_ndims[0]
    else:
        raise ValueError(
            "All coefficients must have a matching number of dimensions")

    if np.isscalar(axes):
        axes = (axes, )
    if axes is None:
        axes = range(ndim)
    else:
        axes = tuple(axes)
    if len(axes) != len(set(axes)):
        raise ValueError("The axes passed to waverecn must be unique.")
    ndim_transform = len(axes)

    for idx, d in enumerate(ds):
        if a is None and not d:
            continue
        # The following if statement handles the case where the approximation
        # coefficient returned at the previous level may exceed the size of the
        # stored detail coefficients by 1 on any given axis.
        if idx > 0:
            a = _match_coeff_dims(a, d)
        d['a' * ndim_transform] = a
        a = idwtn(d, wavelet, mode, axes)

    return a


def _coeffs_wavedec_to_wavedecn(coeffs):
    """Convert wavedec coefficients to the wavedecn format."""
    if len(coeffs) == 0:
        return coeffs
    coeffs = copy(coeffs)
    for n in range(1, len(coeffs)):
        if coeffs[n] is None:
            continue
        if coeffs[n].ndim != 1:
            raise ValueError("expected a 1D coefficient array")
        coeffs[n] = dict(d=coeffs[n])
    return coeffs


def _coeffs_wavedec2_to_wavedecn(coeffs):
    """Convert wavedec2 coefficients to the wavedecn format."""
    if len(coeffs) == 0:
        return coeffs
    coeffs = copy(coeffs)
    for n in range(1, len(coeffs)):
        if not isinstance(coeffs[n], (tuple, list)) or len(coeffs[n]) != 3:
            raise ValueError("expected a 3-tuple of detail coefficients")
        (da, ad, dd) = coeffs[n]
        coeffs[n] = dict(ad=ad, da=da, dd=dd)
    return coeffs


def _determine_coeff_array_shape(coeffs, axes):
    arr_shape = np.asarray(coeffs[0].shape)
    axes = np.asarray(axes)  # axes that were transformed
    ndim_transform = len(axes)
    ncoeffs = coeffs[0].size
    for d in coeffs[1:]:
        arr_shape[axes] += np.asarray(d['d'*ndim_transform].shape)[axes]
        for k, v in d.items():
            ncoeffs += v.size
    arr_shape = tuple(arr_shape.tolist())
    # if the total number of coefficients doesn't equal the size of the array
    # then tight packing is not possible.
    is_tight_packing = (np.prod(arr_shape) == ncoeffs)
    return arr_shape, is_tight_packing


def coeffs_to_array(coeffs, padding=0, axes=None):
    """
    Arrange a wavelet coefficient list from `wavedecn` into a single array.

    Parameters
    ----------

    coeffs: array-like
        dictionary of wavelet coefficients as returned by pywt.wavedecn
    padding : float or None, optional
        If None, raise an error if the coefficients cannot be tightly packed.
    axes : sequence of ints, optional
        Axes over which the DWT that created ``coeffs`` was performed.  The
        default value of ``None`` corresponds to all axes.

    Returns
    -------
    coeff_arr : array-like
        Wavelet transform coefficient array.
    coeff_slices : list
        List of slices corresponding to each coefficient.  As a 2D example,
        `coeff_arr[coeff_slices[1]['dd']]` would extract the first level detail
        coefficients from `coeff_arr`.

    See Also
    --------
    array_to_coeffs : the inverse of coeffs_to_array

    Notes
    -----
    Assume a 2D coefficient dictionary, c, from a two-level transform.

    Then all 2D coefficients will be stacked into a single larger 2D array
    as follows::

        +---------------+---------------+-------------------------------+
        |               |               |                               |
        |     c[0]      |  c[1]['da']   |                               |
        |               |               |                               |
        +---------------+---------------+           c[2]['da']          |
        |               |               |                               |
        | c[1]['ad']    |  c[1]['dd']   |                               |
        |               |               |                               |
        +---------------+---------------+ ------------------------------+
        |                               |                               |
        |                               |                               |
        |                               |                               |
        |          c[2]['ad']           |           c[2]['dd']          |
        |                               |                               |
        |                               |                               |
        |                               |                               |
        +-------------------------------+-------------------------------+

    Examples
    --------
    >>> import pywt
    >>> cam = pywt.data.camera()
    >>> coeffs = pywt.wavedecn(cam, wavelet='db2', level=3)
    >>> arr, coeff_slices = pywt.coeffs_to_array(coeffs)

    """
    if not isinstance(coeffs, list) or len(coeffs) == 0:
        raise ValueError("input must be a list of coefficients from wavedecn")
    if coeffs[0] is None:
        raise ValueError("coeffs_to_array does not support missing "
                         "coefficients.")
    if not isinstance(coeffs[0], np.ndarray):
        raise ValueError("first list element must be a numpy array")
    if len(coeffs) > 1:
        # convert wavedec or wavedec2 format coefficients to waverecn format
        if isinstance(coeffs[1], dict):
            pass
        elif isinstance(coeffs[1], np.ndarray):
            coeffs = _coeffs_wavedec_to_wavedecn(coeffs)
        elif isinstance(coeffs[1], (tuple, list)):
            coeffs = _coeffs_wavedec2_to_wavedecn(coeffs)
        else:
            raise ValueError("invalid coefficient list")
    # initialize with the approximation coefficients.
    a_coeffs = coeffs[0]
    a_shape = a_coeffs.shape
    ndim = a_coeffs.ndim
    if len(coeffs) == 1:
        # only a single approximation coefficient array was found
        return a_coeffs, [[slice(None)] * ndim]

    # Determine the number of dimensions that were transformed via key length
    ndim_transform = len(list(coeffs[1].keys())[0])
    if axes is None:
        if ndim_transform < ndim:
            raise ValueError(
                "coeffs corresponds to a DWT performed over only a subset of "
                "the axes.  In this case, axes must be specified.")
        axes = np.arange(ndim)

    if len(axes) != ndim_transform:
        raise ValueError(
            "The length of axes doesn't match the number of dimensions "
            "transformed.")

    # determine size of output and if tight packing is possible
    arr_shape, is_tight_packing = _determine_coeff_array_shape(coeffs, axes)

    # preallocate output array
    if padding is None:
        if not is_tight_packing:
            raise ValueError("array coefficients cannot be tightly packed")
        coeff_arr = np.empty(arr_shape, dtype=a_coeffs.dtype)
    else:
        coeff_arr = np.full(arr_shape, padding, dtype=a_coeffs.dtype)

    a_slices = [slice(s) for s in a_shape]
    coeff_arr[a_slices] = a_coeffs

    # initialize list of coefficient slices
    coeff_slices = []
    coeff_slices.append(a_slices)

    # loop over the detail cofficients, adding them to coeff_arr
    ds = coeffs[1:]
    for coeff_dict in ds:
        coeff_slices.append({})  # new dictionary for detail coefficients
        if np.any([d is None for d in coeff_dict.values()]):
            raise ValueError("coeffs_to_array does not support missing "
                             "coefficients.")
        d_shape = coeff_dict['d' * ndim_transform].shape
        for key in coeff_dict.keys():
            d = coeff_dict[key]
            slice_array = [slice(None), ] * ndim
            for i, let in enumerate(key):
                ax_i = axes[i]  # axis corresponding to this transform index
                if let == 'a':
                    slice_array[ax_i] = slice(d.shape[ax_i])
                elif let == 'd':
                    slice_array[ax_i] = slice(a_shape[ax_i],
                                              a_shape[ax_i] + d.shape[ax_i])
                else:
                    raise ValueError("unexpected letter: {}".format(let))
            coeff_arr[slice_array] = d
            coeff_slices[-1][key] = slice_array
        a_shape = [a_shape[n] + d_shape[n] for n in range(ndim)]
    return coeff_arr, coeff_slices


def array_to_coeffs(arr, coeff_slices, output_format='wavedecn'):
    """
    Convert a combined array of coefficients back to a list compatible with
    `waverecn`.

    Parameters
    ----------

    arr: array-like
        An array containing all wavelet coefficients.  This should have been
        generated via `coeffs_to_array`.
    coeff_slices : list of tuples
        List of slices corresponding to each coefficient as obtained from
        `array_to_coeffs`.
    output_format : {'wavedec', 'wavedec2', 'wavedecn'}
        Make the form of the coefficients compatible with this type of
        multilevel transform.

    Returns
    -------
    coeffs: array-like
        Wavelet transform coefficient array.

    See Also
    --------
    coeffs_to_array : the inverse of array_to_coeffs

    Notes
    -----
    A single large array containing all coefficients will have subsets stored,
    into a `waverecn` list, c, as indicated below::

        +---------------+---------------+-------------------------------+
        |               |               |                               |
        |     c[0]      |  c[1]['da']   |                               |
        |               |               |                               |
        +---------------+---------------+           c[2]['da']          |
        |               |               |                               |
        | c[1]['ad']    |  c[1]['dd']   |                               |
        |               |               |                               |
        +---------------+---------------+ ------------------------------+
        |                               |                               |
        |                               |                               |
        |                               |                               |
        |          c[2]['ad']           |           c[2]['dd']          |
        |                               |                               |
        |                               |                               |
        |                               |                               |
        +-------------------------------+-------------------------------+

    Examples
    --------
    >>> import pywt
    >>> from numpy.testing import assert_array_almost_equal
    >>> cam = pywt.data.camera()
    >>> coeffs = pywt.wavedecn(cam, wavelet='db2', level=3)
    >>> arr, coeff_slices = pywt.coeffs_to_array(coeffs)
    >>> coeffs_from_arr = pywt.array_to_coeffs(arr, coeff_slices)
    >>> cam_recon = pywt.waverecn(coeffs_from_arr, wavelet='db2')
    >>> assert_array_almost_equal(cam, cam_recon)

    """
    arr = np.asarray(arr)
    coeffs = []
    if len(coeff_slices) == 0:
        raise ValueError("empty list of coefficient slices")
    else:
        coeffs.append(arr[coeff_slices[0]])

    # difference coefficients at each level
    for n in range(1, len(coeff_slices)):
        if output_format == 'wavedec':
            d = arr[coeff_slices[n]['d']]
        elif output_format == 'wavedec2':
            d = (arr[coeff_slices[n]['da']],
                 arr[coeff_slices[n]['ad']],
                 arr[coeff_slices[n]['dd']])
        elif output_format == 'wavedecn':
            d = {}
            for k, v in coeff_slices[n].items():
                d[k] = arr[v]
        else:
            raise ValueError(
                "Unrecognized output format: {}".format(output_format))
        coeffs.append(d)
    return coeffs
