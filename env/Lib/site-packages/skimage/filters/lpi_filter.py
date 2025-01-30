"""
:author: Stefan van der Walt, 2008
:license: modified BSD
"""

import numpy as np
import scipy.fft as fft

from .._shared.utils import _supported_float_type, check_nD


def _min_limit(x, val=np.finfo(float).eps):
    mask = np.abs(x) < val
    x[mask] = np.sign(x[mask]) * val


def _center(x, oshape):
    """Return an array of shape ``oshape`` from the center of array ``x``."""
    start = (np.array(x.shape) - np.array(oshape)) // 2
    out = x[tuple(slice(s, s + n) for s, n in zip(start, oshape))]
    return out


def _pad(data, shape):
    """Pad the data to the given shape with zeros.

    Parameters
    ----------
    data : 2-d ndarray
        Input data
    shape : (2,) tuple

    """
    out = np.zeros(shape, dtype=data.dtype)
    out[tuple(slice(0, n) for n in data.shape)] = data
    return out


class LPIFilter2D:
    """Linear Position-Invariant Filter (2-dimensional)"""

    def __init__(self, impulse_response, **filter_params):
        """
        Parameters
        ----------
        impulse_response : callable `f(r, c, **filter_params)`
            Function that yields the impulse response.  ``r`` and ``c`` are
            1-dimensional vectors that represent row and column positions, in
            other words coordinates are (r[0],c[0]),(r[0],c[1]) etc.
            `**filter_params` are passed through.

            In other words, ``impulse_response`` would be called like this:

            >>> def impulse_response(r, c, **filter_params):
            ...     pass
            >>>
            >>> r = [0,0,0,1,1,1,2,2,2]
            >>> c = [0,1,2,0,1,2,0,1,2]
            >>> filter_params = {'kw1': 1, 'kw2': 2, 'kw3': 3}
            >>> impulse_response(r, c, **filter_params)


        Examples
        --------
        Gaussian filter without normalization of coefficients:

        >>> def filt_func(r, c, sigma=1):
        ...     return np.exp(-(r**2 + c**2)/(2 * sigma**2))
        >>> filter = LPIFilter2D(filt_func)

        """
        if not callable(impulse_response):
            raise ValueError("Impulse response must be a callable.")

        self.impulse_response = impulse_response
        self.filter_params = filter_params
        self._cache = None

    def _prepare(self, data):
        """Calculate filter and data FFT in preparation for filtering."""
        dshape = np.array(data.shape)
        even_offset = (dshape % 2 == 0).astype(int)
        dshape += even_offset  # all filter dimensions must be uneven
        oshape = np.array(data.shape) * 2 - 1

        float_dtype = _supported_float_type(data.dtype)
        data = data.astype(float_dtype, copy=False)

        if self._cache is None or np.any(self._cache.shape != oshape):
            coords = np.mgrid[
                [
                    slice(0 + offset, float(n + offset))
                    for (n, offset) in zip(dshape, even_offset)
                ]
            ]
            # this steps over two sets of coordinates,
            # not over the coordinates individually
            for k, coord in enumerate(coords):
                coord -= (dshape[k] - 1) / 2.0
            coords = coords.reshape(2, -1).T  # coordinate pairs (r,c)
            coords = coords.astype(float_dtype, copy=False)

            f = self.impulse_response(
                coords[:, 0], coords[:, 1], **self.filter_params
            ).reshape(dshape)

            f = _pad(f, oshape)
            F = fft.fftn(f)
            self._cache = F
        else:
            F = self._cache

        data = _pad(data, oshape)
        G = fft.fftn(data)

        return F, G

    def __call__(self, data):
        """Apply the filter to the given data.

        Parameters
        ----------
        data : (M, N) ndarray

        """
        check_nD(data, 2, 'data')
        F, G = self._prepare(data)
        out = fft.ifftn(F * G)
        out = np.abs(_center(out, data.shape))
        return out


def filter_forward(
    data, impulse_response=None, filter_params=None, predefined_filter=None
):
    """Apply the given filter to data.

    Parameters
    ----------
    data : (M, N) ndarray
        Input data.
    impulse_response : callable `f(r, c, **filter_params)`
        Impulse response of the filter.  See LPIFilter2D.__init__.
    filter_params : dict, optional
        Additional keyword parameters to the impulse_response function.

    Other Parameters
    ----------------
    predefined_filter : LPIFilter2D
        If you need to apply the same filter multiple times over different
        images, construct the LPIFilter2D and specify it here.

    Examples
    --------

    Gaussian filter without normalization:

    >>> def filt_func(r, c, sigma=1):
    ...     return np.exp(-(r**2 + c**2)/(2 * sigma**2))
    >>>
    >>> from skimage import data
    >>> filtered = filter_forward(data.coins(), filt_func)

    """
    if filter_params is None:
        filter_params = {}
    check_nD(data, 2, 'data')
    if predefined_filter is None:
        predefined_filter = LPIFilter2D(impulse_response, **filter_params)
    return predefined_filter(data)


def filter_inverse(
    data, impulse_response=None, filter_params=None, max_gain=2, predefined_filter=None
):
    """Apply the filter in reverse to the given data.

    Parameters
    ----------
    data : (M, N) ndarray
        Input data.
    impulse_response : callable `f(r, c, **filter_params)`
        Impulse response of the filter.  See :class:`~.LPIFilter2D`. This is a required
        argument unless a `predifined_filter` is provided.
    filter_params : dict, optional
        Additional keyword parameters to the impulse_response function.
    max_gain : float, optional
        Limit the filter gain.  Often, the filter contains zeros, which would
        cause the inverse filter to have infinite gain.  High gain causes
        amplification of artefacts, so a conservative limit is recommended.

    Other Parameters
    ----------------
    predefined_filter : LPIFilter2D, optional
        If you need to apply the same filter multiple times over different
        images, construct the LPIFilter2D and specify it here.

    """
    if filter_params is None:
        filter_params = {}

    check_nD(data, 2, 'data')
    if predefined_filter is None:
        filt = LPIFilter2D(impulse_response, **filter_params)
    else:
        filt = predefined_filter

    F, G = filt._prepare(data)
    _min_limit(F, val=np.finfo(F.real.dtype).eps)

    F = 1 / F
    mask = np.abs(F) > max_gain
    F[mask] = np.sign(F[mask]) * max_gain

    return _center(np.abs(fft.ifftshift(fft.ifftn(G * F))), data.shape)


def wiener(
    data, impulse_response=None, filter_params=None, K=0.25, predefined_filter=None
):
    """Minimum Mean Square Error (Wiener) inverse filter.

    Parameters
    ----------
    data : (M, N) ndarray
        Input data.
    K : float or (M, N) ndarray
        Ratio between power spectrum of noise and undegraded
        image.
    impulse_response : callable `f(r, c, **filter_params)`
        Impulse response of the filter.  See LPIFilter2D.__init__.
    filter_params : dict, optional
        Additional keyword parameters to the impulse_response function.

    Other Parameters
    ----------------
    predefined_filter : LPIFilter2D
        If you need to apply the same filter multiple times over different
        images, construct the LPIFilter2D and specify it here.

    """
    if filter_params is None:
        filter_params = {}

    check_nD(data, 2, 'data')

    if not isinstance(K, float):
        check_nD(K, 2, 'K')

    if predefined_filter is None:
        filt = LPIFilter2D(impulse_response, **filter_params)
    else:
        filt = predefined_filter

    F, G = filt._prepare(data)
    _min_limit(F, val=np.finfo(F.real.dtype).eps)

    H_mag_sqr = np.abs(F) ** 2
    F = 1 / F * H_mag_sqr / (H_mag_sqr + K)

    return _center(np.abs(fft.ifftshift(fft.ifftn(G * F))), data.shape)
