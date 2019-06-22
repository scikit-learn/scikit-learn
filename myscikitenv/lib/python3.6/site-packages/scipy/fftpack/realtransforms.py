"""
Real spectrum transforms (DCT, DST, MDCT)
"""
from __future__ import division, print_function, absolute_import


__all__ = ['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn', 'idstn']

import numpy as np
from scipy.fftpack import _fftpack
from scipy.fftpack.basic import _datacopied, _fix_shape, _asfarray
from scipy.fftpack.helper import _init_nd_shape_and_axes

import atexit
atexit.register(_fftpack.destroy_ddct1_cache)
atexit.register(_fftpack.destroy_ddct2_cache)
atexit.register(_fftpack.destroy_ddct4_cache)
atexit.register(_fftpack.destroy_dct1_cache)
atexit.register(_fftpack.destroy_dct2_cache)
atexit.register(_fftpack.destroy_dct4_cache)

atexit.register(_fftpack.destroy_ddst1_cache)
atexit.register(_fftpack.destroy_ddst2_cache)
atexit.register(_fftpack.destroy_dst1_cache)
atexit.register(_fftpack.destroy_dst2_cache)


def dctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return multidimensional Discrete Cosine Transform along the specified axes.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DCT (see Notes). Default type is 2.
    shape : int or array_like of ints or None, optional
        The shape of the result.  If both `shape` and `axes` (see below) are
        None, `shape` is ``x.shape``; if `shape` is None but `axes` is
        not None, then `shape` is ``scipy.take(x.shape, axes, axis=0)``.
        If ``shape[i] > x.shape[i]``, the i-th dimension is padded with zeros.
        If ``shape[i] < x.shape[i]``, the i-th dimension is truncated to
        length ``shape[i]``.
        If any element of `shape` is -1, the size of the corresponding
        dimension of `x` is used.
    axes : int or array_like of ints or None, optional
        Axes along which the DCT is computed.
        The default is over all axes.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    y : ndarray of real
        The transformed input array.

    See Also
    --------
    idctn : Inverse multidimensional DCT

    Notes
    -----
    For full details of the DCT types and normalization modes, as well as
    references, see `dct`.

    Examples
    --------
    >>> from scipy.fftpack import dctn, idctn
    >>> y = np.random.randn(16, 16)
    >>> np.allclose(y, idctn(dctn(y, norm='ortho'), norm='ortho'))
    True

    """
    x = np.asanyarray(x)
    shape, axes = _init_nd_shape_and_axes(x, shape, axes)
    for n, ax in zip(shape, axes):
        x = dct(x, type=type, n=n, axis=ax, norm=norm, overwrite_x=overwrite_x)
    return x


def idctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return multidimensional Discrete Cosine Transform along the specified axes.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DCT (see Notes). Default type is 2.
    shape : int or array_like of ints or None, optional
        The shape of the result.  If both `shape` and `axes` (see below) are
        None, `shape` is ``x.shape``; if `shape` is None but `axes` is
        not None, then `shape` is ``scipy.take(x.shape, axes, axis=0)``.
        If ``shape[i] > x.shape[i]``, the i-th dimension is padded with zeros.
        If ``shape[i] < x.shape[i]``, the i-th dimension is truncated to
        length ``shape[i]``.
        If any element of `shape` is -1, the size of the corresponding
        dimension of `x` is used.
    axes : int or array_like of ints or None, optional
        Axes along which the IDCT is computed.
        The default is over all axes.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    y : ndarray of real
        The transformed input array.

    See Also
    --------
    dctn : multidimensional DCT

    Notes
    -----
    For full details of the IDCT types and normalization modes, as well as
    references, see `idct`.

    Examples
    --------
    >>> from scipy.fftpack import dctn, idctn
    >>> y = np.random.randn(16, 16)
    >>> np.allclose(y, idctn(dctn(y, norm='ortho'), norm='ortho'))
    True

    """
    x = np.asanyarray(x)
    shape, axes = _init_nd_shape_and_axes(x, shape, axes)
    for n, ax in zip(shape, axes):
        x = idct(x, type=type, n=n, axis=ax, norm=norm,
                 overwrite_x=overwrite_x)
    return x


def dstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return multidimensional Discrete Sine Transform along the specified axes.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DST (see Notes). Default type is 2.
    shape : int or array_like of ints or None, optional
        The shape of the result.  If both `shape` and `axes` (see below) are
        None, `shape` is ``x.shape``; if `shape` is None but `axes` is
        not None, then `shape` is ``scipy.take(x.shape, axes, axis=0)``.
        If ``shape[i] > x.shape[i]``, the i-th dimension is padded with zeros.
        If ``shape[i] < x.shape[i]``, the i-th dimension is truncated to
        length ``shape[i]``.
        If any element of `shape` is -1, the size of the corresponding
        dimension of `x` is used.
    axes : int or array_like of ints or None, optional
        Axes along which the DCT is computed.
        The default is over all axes.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    y : ndarray of real
        The transformed input array.

    See Also
    --------
    idstn : Inverse multidimensional DST

    Notes
    -----
    For full details of the DST types and normalization modes, as well as
    references, see `dst`.

    Examples
    --------
    >>> from scipy.fftpack import dstn, idstn
    >>> y = np.random.randn(16, 16)
    >>> np.allclose(y, idstn(dstn(y, norm='ortho'), norm='ortho'))
    True

    """
    x = np.asanyarray(x)
    shape, axes = _init_nd_shape_and_axes(x, shape, axes)
    for n, ax in zip(shape, axes):
        x = dst(x, type=type, n=n, axis=ax, norm=norm, overwrite_x=overwrite_x)
    return x


def idstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return multidimensional Discrete Sine Transform along the specified axes.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DST (see Notes). Default type is 2.
    shape : int or array_like of ints or None, optional
        The shape of the result.  If both `shape` and `axes` (see below) are
        None, `shape` is ``x.shape``; if `shape` is None but `axes` is
        not None, then `shape` is ``scipy.take(x.shape, axes, axis=0)``.
        If ``shape[i] > x.shape[i]``, the i-th dimension is padded with zeros.
        If ``shape[i] < x.shape[i]``, the i-th dimension is truncated to
        length ``shape[i]``.
        If any element of `shape` is -1, the size of the corresponding
        dimension of `x` is used.
    axes : int or array_like of ints or None, optional
        Axes along which the IDST is computed.
        The default is over all axes.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    y : ndarray of real
        The transformed input array.

    See Also
    --------
    dstn : multidimensional DST

    Notes
    -----
    For full details of the IDST types and normalization modes, as well as
    references, see `idst`.

    Examples
    --------
    >>> from scipy.fftpack import dstn, idstn
    >>> y = np.random.randn(16, 16)
    >>> np.allclose(y, idstn(dstn(y, norm='ortho'), norm='ortho'))
    True

    """
    x = np.asanyarray(x)
    shape, axes = _init_nd_shape_and_axes(x, shape, axes)
    for n, ax in zip(shape, axes):
        x = idst(x, type=type, n=n, axis=ax, norm=norm,
                 overwrite_x=overwrite_x)
    return x


def dct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return the Discrete Cosine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DCT (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.  If ``n < x.shape[axis]``, `x` is
        truncated.  If ``n > x.shape[axis]``, `x` is zero-padded. The
        default results in ``n = x.shape[axis]``.
    axis : int, optional
        Axis along which the dct is computed; the default is over the
        last axis (i.e., ``axis=-1``).
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    y : ndarray of real
        The transformed input array.

    See Also
    --------
    idct : Inverse DCT

    Notes
    -----
    For a single dimension array ``x``, ``dct(x, norm='ortho')`` is equal to
    MATLAB ``dct(x)``.

    There are theoretically 8 types of the DCT, only the first 4 types are
    implemented in scipy. 'The' DCT generally refers to DCT type 2, and 'the'
    Inverse DCT generally refers to DCT type 3.

    **Type I**

    There are several definitions of the DCT-I; we use the following
    (for ``norm=None``)::

                                         N-2
      y[k] = x[0] + (-1)**k x[N-1] + 2 * sum x[n]*cos(pi*k*n/(N-1))
                                         n=1

    If ``norm='ortho'``, ``x[0]`` and ``x[N-1]`` are multiplied by
    a scaling factor of ``sqrt(2)``, and ``y[k]`` is multiplied by a
    scaling factor `f`::

      f = 0.5*sqrt(1/(N-1)) if k = 0 or N-1,
      f = 0.5*sqrt(2/(N-1)) otherwise.

    .. versionadded:: 1.2.0
       Orthonormalization in DCT-I.

    .. note::
       The DCT-I is only supported for input size > 1.

    **Type II**

    There are several definitions of the DCT-II; we use the following
    (for ``norm=None``)::


                N-1
      y[k] = 2* sum x[n]*cos(pi*k*(2n+1)/(2*N)), 0 <= k < N.
                n=0

    If ``norm='ortho'``, ``y[k]`` is multiplied by a scaling factor `f`::

      f = sqrt(1/(4*N)) if k = 0,
      f = sqrt(1/(2*N)) otherwise.

    Which makes the corresponding matrix of coefficients orthonormal
    (``OO' = Id``).

    **Type III**

    There are several definitions, we use the following
    (for ``norm=None``)::

                        N-1
      y[k] = x[0] + 2 * sum x[n]*cos(pi*(k+0.5)*n/N), 0 <= k < N.
                        n=1

    or, for ``norm='ortho'`` and 0 <= k < N::

                                          N-1
      y[k] = x[0] / sqrt(N) + sqrt(2/N) * sum x[n]*cos(pi*(k+0.5)*n/N)
                                          n=1

    The (unnormalized) DCT-III is the inverse of the (unnormalized) DCT-II, up
    to a factor `2N`. The orthonormalized DCT-III is exactly the inverse of
    the orthonormalized DCT-II.

    **Type IV**

    There are several definitions of the DCT-IV; we use the following
    (for ``norm=None``)::


                N-1
      y[k] = 2* sum x[n]*cos(pi*(2k+1)*(2n+1)/(4*N)), 0 <= k < N.
                n=0

    If ``norm='ortho'``, ``y[k]`` is multiplied by a scaling factor `f`::

      f = 0.5*sqrt(2/N)

    .. versionadded:: 1.2.0
       Support for DCT-IV.

    References
    ----------
    .. [1] 'A Fast Cosine Transform in One and Two Dimensions', by J.
           Makhoul, `IEEE Transactions on acoustics, speech and signal
           processing` vol. 28(1), pp. 27-34,
           :doi:`10.1109/TASSP.1980.1163351` (1980).
    .. [2] Wikipedia, "Discrete cosine transform",
           https://en.wikipedia.org/wiki/Discrete_cosine_transform

    Examples
    --------
    The Type 1 DCT is equivalent to the FFT (though faster) for real,
    even-symmetrical inputs.  The output is also real and even-symmetrical.
    Half of the FFT input is used to generate half of the FFT output:

    >>> from scipy.fftpack import fft, dct
    >>> fft(np.array([4., 3., 5., 10., 5., 3.])).real
    array([ 30.,  -8.,   6.,  -2.,   6.,  -8.])
    >>> dct(np.array([4., 3., 5., 10.]), 1)
    array([ 30.,  -8.,   6.,  -2.])

    """
    return _dct(x, type, n, axis, normalize=norm, overwrite_x=overwrite_x)


def idct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return the Inverse Discrete Cosine Transform of an arbitrary type sequence.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DCT (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.  If ``n < x.shape[axis]``, `x` is
        truncated.  If ``n > x.shape[axis]``, `x` is zero-padded. The
        default results in ``n = x.shape[axis]``.
    axis : int, optional
        Axis along which the idct is computed; the default is over the
        last axis (i.e., ``axis=-1``).
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    idct : ndarray of real
        The transformed input array.

    See Also
    --------
    dct : Forward DCT

    Notes
    -----
    For a single dimension array `x`, ``idct(x, norm='ortho')`` is equal to
    MATLAB ``idct(x)``.

    'The' IDCT is the IDCT of type 2, which is the same as DCT of type 3.

    IDCT of type 1 is the DCT of type 1, IDCT of type 2 is the DCT of type
    3, and IDCT of type 3 is the DCT of type 2. IDCT of type 4 is the DCT
    of type 4. For the definition of these types, see `dct`.

    Examples
    --------
    The Type 1 DCT is equivalent to the DFT for real, even-symmetrical
    inputs.  The output is also real and even-symmetrical.  Half of the IFFT
    input is used to generate half of the IFFT output:

    >>> from scipy.fftpack import ifft, idct
    >>> ifft(np.array([ 30.,  -8.,   6.,  -2.,   6.,  -8.])).real
    array([  4.,   3.,   5.,  10.,   5.,   3.])
    >>> idct(np.array([ 30.,  -8.,   6.,  -2.]), 1) / 6
    array([  4.,   3.,   5.,  10.])

    """
    # Inverse/forward type table
    _TP = {1:1, 2:3, 3:2, 4:4}
    return _dct(x, _TP[type], n, axis, normalize=norm, overwrite_x=overwrite_x)


def _get_dct_fun(type, dtype):
    try:
        name = {'float64':'ddct%d', 'float32':'dct%d'}[dtype.name]
    except KeyError:
        raise ValueError("dtype %s not supported" % dtype)
    try:
        f = getattr(_fftpack, name % type)
    except AttributeError as e:
        raise ValueError(str(e) + ". Type %d not understood" % type)
    return f


def _get_norm_mode(normalize):
    try:
        nm = {None:0, 'ortho':1}[normalize]
    except KeyError:
        raise ValueError("Unknown normalize mode %s" % normalize)
    return nm


def __fix_shape(x, n, axis, dct_or_dst):
    tmp = _asfarray(x)
    copy_made = _datacopied(tmp, x)
    if n is None:
        n = tmp.shape[axis]
    elif n != tmp.shape[axis]:
        tmp, copy_made2 = _fix_shape(tmp, n, axis)
        copy_made = copy_made or copy_made2
    if n < 1:
        raise ValueError("Invalid number of %s data points "
                         "(%d) specified." % (dct_or_dst, n))
    return tmp, n, copy_made


def _raw_dct(x0, type, n, axis, nm, overwrite_x):
    f = _get_dct_fun(type, x0.dtype)
    return _eval_fun(f, x0, n, axis, nm, overwrite_x)


def _raw_dst(x0, type, n, axis, nm, overwrite_x):
    f = _get_dst_fun(type, x0.dtype)
    return _eval_fun(f, x0, n, axis, nm, overwrite_x)


def _eval_fun(f, tmp, n, axis, nm, overwrite_x):
    if axis == -1 or axis == len(tmp.shape) - 1:
        return f(tmp, n, nm, overwrite_x)

    tmp = np.swapaxes(tmp, axis, -1)
    tmp = f(tmp, n, nm, overwrite_x)
    return np.swapaxes(tmp, axis, -1)


def _dct(x, type, n=None, axis=-1, overwrite_x=False, normalize=None):
    """
    Return Discrete Cosine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array_like
        input array.
    n : int, optional
        Length of the transform.  If ``n < x.shape[axis]``, `x` is
        truncated.  If ``n > x.shape[axis]``, `x` is zero-padded. The
        default results in ``n = x.shape[axis]``.
    axis : int, optional
        Axis along which the dct is computed; the default is over the
        last axis (i.e., ``axis=-1``).
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    z : ndarray

    """
    x0, n, copy_made = __fix_shape(x, n, axis, 'DCT')
    if type == 1 and n < 2:
        raise ValueError("DCT-I is not defined for size < 2")
    overwrite_x = overwrite_x or copy_made
    nm = _get_norm_mode(normalize)
    if np.iscomplexobj(x0):
        return (_raw_dct(x0.real, type, n, axis, nm, overwrite_x) + 1j *
                _raw_dct(x0.imag, type, n, axis, nm, overwrite_x))
    else:
        return _raw_dct(x0, type, n, axis, nm, overwrite_x)


def dst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return the Discrete Sine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DST (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.  If ``n < x.shape[axis]``, `x` is
        truncated.  If ``n > x.shape[axis]``, `x` is zero-padded. The
        default results in ``n = x.shape[axis]``.
    axis : int, optional
        Axis along which the dst is computed; the default is over the
        last axis (i.e., ``axis=-1``).
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    dst : ndarray of reals
        The transformed input array.

    See Also
    --------
    idst : Inverse DST

    Notes
    -----
    For a single dimension array ``x``.

    There are theoretically 8 types of the DST for different combinations of
    even/odd boundary conditions and boundary off sets [1]_, only the first
    4 types are implemented in scipy.

    **Type I**

    There are several definitions of the DST-I; we use the following
    for ``norm=None``.  DST-I assumes the input is odd around n=-1 and n=N. ::

                 N-1
      y[k] = 2 * sum x[n]*sin(pi*(k+1)*(n+1)/(N+1))
                 n=0

    Note that the DST-I is only supported for input size > 1
    The (unnormalized) DST-I is its own inverse, up to a factor `2(N+1)`.
    The orthonormalized DST-I is exactly its own inverse.

    **Type II**

    There are several definitions of the DST-II; we use the following
    for ``norm=None``.  DST-II assumes the input is odd around n=-1/2 and
    n=N-1/2; the output is odd around k=-1 and even around k=N-1 ::

                N-1
      y[k] = 2* sum x[n]*sin(pi*(k+1)*(n+0.5)/N), 0 <= k < N.
                n=0

    if ``norm='ortho'``, ``y[k]`` is multiplied by a scaling factor `f` ::

        f = sqrt(1/(4*N)) if k == 0
        f = sqrt(1/(2*N)) otherwise.

    **Type III**

    There are several definitions of the DST-III, we use the following
    (for ``norm=None``).  DST-III assumes the input is odd around n=-1
    and even around n=N-1 ::

                                 N-2
      y[k] = x[N-1]*(-1)**k + 2* sum x[n]*sin(pi*(k+0.5)*(n+1)/N), 0 <= k < N.
                                 n=0

    The (unnormalized) DST-III is the inverse of the (unnormalized) DST-II, up
    to a factor `2N`.  The orthonormalized DST-III is exactly the inverse of
    the orthonormalized DST-II.

    .. versionadded:: 0.11.0

    **Type IV**

    There are several definitions of the DST-IV, we use the following
    (for ``norm=None``).  DST-IV assumes the input is odd around n=-0.5
    and even around n=N-0.5 ::

                N-1
      y[k] = 2* sum x[n]*sin(pi*(k+0.5)*(n+0.5)/N), 0 <= k < N.
                n=0

    The (unnormalized) DST-IV is its own inverse, up
    to a factor `2N`.  The orthonormalized DST-IV is exactly its own inverse.

    .. versionadded:: 1.2.0
       Support for DST-IV.

    References
    ----------
    .. [1] Wikipedia, "Discrete sine transform",
           https://en.wikipedia.org/wiki/Discrete_sine_transform

    """
    return _dst(x, type, n, axis, normalize=norm, overwrite_x=overwrite_x)


def idst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return the Inverse Discrete Sine Transform of an arbitrary type sequence.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3, 4}, optional
        Type of the DST (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.  If ``n < x.shape[axis]``, `x` is
        truncated.  If ``n > x.shape[axis]``, `x` is zero-padded. The
        default results in ``n = x.shape[axis]``.
    axis : int, optional
        Axis along which the idst is computed; the default is over the
        last axis (i.e., ``axis=-1``).
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True, the contents of `x` can be destroyed; the default is False.

    Returns
    -------
    idst : ndarray of real
        The transformed input array.

    See Also
    --------
    dst : Forward DST

    Notes
    -----
    'The' IDST is the IDST of type 2, which is the same as DST of type 3.

    IDST of type 1 is the DST of type 1, IDST of type 2 is the DST of type
    3, and IDST of type 3 is the DST of type 2. For the definition of these
    types, see `dst`.

    .. versionadded:: 0.11.0

    """
    # Inverse/forward type table
    _TP = {1:1, 2:3, 3:2, 4:4}
    return _dst(x, _TP[type], n, axis, normalize=norm, overwrite_x=overwrite_x)


def _get_dst_fun(type, dtype):
    try:
        name = {'float64':'ddst%d', 'float32':'dst%d'}[dtype.name]
    except KeyError:
        raise ValueError("dtype %s not supported" % dtype)
    try:
        f = getattr(_fftpack, name % type)
    except AttributeError as e:
        raise ValueError(str(e) + ". Type %d not understood" % type)
    return f


def _dst(x, type, n=None, axis=-1, overwrite_x=False, normalize=None):
    """
    Return Discrete Sine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array_like
        input array.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis along which the dst is computed. (default=-1)
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    z : real ndarray

    """
    x0, n, copy_made = __fix_shape(x, n, axis, 'DST')
    if type == 1 and n < 2:
        raise ValueError("DST-I is not defined for size < 2")
    overwrite_x = overwrite_x or copy_made
    nm = _get_norm_mode(normalize)
    if np.iscomplexobj(x0):
        return (_raw_dst(x0.real, type, n, axis, nm, overwrite_x) + 1j *
                _raw_dst(x0.imag, type, n, axis, nm, overwrite_x))
    else:
        return _raw_dst(x0, type, n, axis, nm, overwrite_x)
