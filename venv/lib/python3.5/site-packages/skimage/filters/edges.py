"""

Sobel and Prewitt filters originally part of CellProfiler, code licensed under
both GPL and BSD licenses.
Website: http://www.cellprofiler.org
Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2011 Broad Institute
All rights reserved.
Original author: Lee Kamentsky

"""
import numpy as np
from .. import img_as_float
from .._shared.utils import assert_nD
from scipy.ndimage import convolve, binary_erosion, generate_binary_structure

from ..restoration.uft import laplacian

EROSION_SELEM = generate_binary_structure(2, 2)

HSOBEL_WEIGHTS = np.array([[ 1, 2, 1],
                           [ 0, 0, 0],
                           [-1,-2,-1]]) / 4.0
VSOBEL_WEIGHTS = HSOBEL_WEIGHTS.T

HSCHARR_WEIGHTS = np.array([[ 3,  10,  3],
                            [ 0,   0,  0],
                            [-3, -10, -3]]) / 16.0
VSCHARR_WEIGHTS = HSCHARR_WEIGHTS.T

HPREWITT_WEIGHTS = np.array([[ 1, 1, 1],
                             [ 0, 0, 0],
                             [-1,-1,-1]]) / 3.0
VPREWITT_WEIGHTS = HPREWITT_WEIGHTS.T

ROBERTS_PD_WEIGHTS = np.array([[1, 0],
                               [0, -1]], dtype=np.double)
ROBERTS_ND_WEIGHTS = np.array([[0, 1],
                               [-1, 0]], dtype=np.double)


def _mask_filter_result(result, mask):
    """Return result after masking.

    Input masks are eroded so that mask areas in the original image don't
    affect values in the result.
    """
    if mask is None:
        result[0, :] = 0
        result[-1, :] = 0
        result[:, 0] = 0
        result[:, -1] = 0
        return result
    else:
        mask = binary_erosion(mask, EROSION_SELEM, border_value=0)
        return result * mask


def sobel(image, mask=None):
    """Find the edge magnitude using the Sobel transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Sobel edge map.

    See also
    --------
    scharr, prewitt, roberts, feature.canny

    Notes
    -----
    Take the square root of the sum of the squares of the horizontal and
    vertical Sobels to get a magnitude that's somewhat insensitive to
    direction.

    The 3x3 convolution kernel used in the horizontal and vertical Sobels is
    an approximation of the gradient of the image (with some slight blurring
    since 9 pixels are used to compute the gradient at a given pixel). As an
    approximation of the gradient, the Sobel operator is not completely
    rotation-invariant. The Scharr operator should be used for a better
    rotation invariance.

    Note that ``scipy.ndimage.sobel`` returns a directional Sobel which
    has to be further processed to perform edge detection.

    Examples
    --------
    >>> from skimage import data
    >>> camera = data.camera()
    >>> from skimage import filters
    >>> edges = filters.sobel(camera)
    """
    assert_nD(image, 2)
    out = np.sqrt(sobel_h(image, mask)**2 + sobel_v(image, mask)**2)
    out /= np.sqrt(2)
    return out


def sobel_h(image, mask=None):
    """Find the horizontal edges of an image using the Sobel transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Sobel edge map.

    Notes
    -----
    We use the following kernel::

      1   2   1
      0   0   0
     -1  -2  -1

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, HSOBEL_WEIGHTS)
    return _mask_filter_result(result, mask)


def sobel_v(image, mask=None):
    """Find the vertical edges of an image using the Sobel transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Sobel edge map.

    Notes
    -----
    We use the following kernel::

      1   0  -1
      2   0  -2
      1   0  -1

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, VSOBEL_WEIGHTS)
    return _mask_filter_result(result, mask)


def scharr(image, mask=None):
    """Find the edge magnitude using the Scharr transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Scharr edge map.

    See also
    --------
    sobel, prewitt, canny

    Notes
    -----
    Take the square root of the sum of the squares of the horizontal and
    vertical Scharrs to get a magnitude that is somewhat insensitive to
    direction. The Scharr operator has a better rotation invariance than
    other edge filters such as the Sobel or the Prewitt operators.

    References
    ----------
    .. [1] D. Kroon, 2009, Short Paper University Twente, Numerical
           Optimization of Kernel Based Image Derivatives.

    .. [2] http://en.wikipedia.org/wiki/Sobel_operator#Alternative_operators

    Examples
    --------
    >>> from skimage import data
    >>> camera = data.camera()
    >>> from skimage import filters
    >>> edges = filters.scharr(camera)
    """
    out = np.sqrt(scharr_h(image, mask)**2 + scharr_v(image, mask)**2)
    out /= np.sqrt(2)
    return out


def scharr_h(image, mask=None):
    """Find the horizontal edges of an image using the Scharr transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Scharr edge map.

    Notes
    -----
    We use the following kernel::

      3   10   3
      0    0   0
     -3  -10  -3

    References
    ----------
    .. [1] D. Kroon, 2009, Short Paper University Twente, Numerical
           Optimization of Kernel Based Image Derivatives.

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, HSCHARR_WEIGHTS)
    return _mask_filter_result(result, mask)


def scharr_v(image, mask=None):
    """Find the vertical edges of an image using the Scharr transform.

    Parameters
    ----------
    image : 2-D array
        Image to process
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Scharr edge map.

    Notes
    -----
    We use the following kernel::

       3   0   -3
      10   0  -10
       3   0   -3

    References
    ----------
    .. [1] D. Kroon, 2009, Short Paper University Twente, Numerical
           Optimization of Kernel Based Image Derivatives.

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, VSCHARR_WEIGHTS)
    return _mask_filter_result(result, mask)


def prewitt(image, mask=None):
    """Find the edge magnitude using the Prewitt transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Prewitt edge map.

    See also
    --------
    sobel, scharr

    Notes
    -----
    Return the square root of the sum of squares of the horizontal
    and vertical Prewitt transforms. The edge magnitude depends slightly
    on edge directions, since the approximation of the gradient operator by
    the Prewitt operator is not completely rotation invariant. For a better
    rotation invariance, the Scharr operator should be used. The Sobel operator
    has a better rotation invariance than the Prewitt operator, but a worse
    rotation invariance than the Scharr operator.

    Examples
    --------
    >>> from skimage import data
    >>> camera = data.camera()
    >>> from skimage import filters
    >>> edges = filters.prewitt(camera)
    """
    assert_nD(image, 2)
    out = np.sqrt(prewitt_h(image, mask)**2 + prewitt_v(image, mask)**2)
    out /= np.sqrt(2)
    return out


def prewitt_h(image, mask=None):
    """Find the horizontal edges of an image using the Prewitt transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Prewitt edge map.

    Notes
    -----
    We use the following kernel::

      1   1   1
      0   0   0
     -1  -1  -1

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, HPREWITT_WEIGHTS)
    return _mask_filter_result(result, mask)


def prewitt_v(image, mask=None):
    """Find the vertical edges of an image using the Prewitt transform.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Prewitt edge map.

    Notes
    -----
    We use the following kernel::

      1   0  -1
      1   0  -1
      1   0  -1

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, VPREWITT_WEIGHTS)
    return _mask_filter_result(result, mask)


def roberts(image, mask=None):
    """Find the edge magnitude using Roberts' cross operator.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Roberts' Cross edge map.

    See also
    --------
    sobel, scharr, prewitt, feature.canny

    Examples
    --------
    >>> from skimage import data
    >>> camera = data.camera()
    >>> from skimage import filters
    >>> edges = filters.roberts(camera)

    """
    assert_nD(image, 2)
    out = np.sqrt(roberts_pos_diag(image, mask)**2 +
                  roberts_neg_diag(image, mask)**2)
    out /= np.sqrt(2)
    return out


def roberts_pos_diag(image, mask=None):
    """Find the cross edges of an image using Roberts' cross operator.

    The kernel is applied to the input image to produce separate measurements
    of the gradient component one orientation.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Robert's edge map.

    Notes
    -----
    We use the following kernel::

      1   0
      0  -1

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, ROBERTS_PD_WEIGHTS)
    return _mask_filter_result(result, mask)


def roberts_neg_diag(image, mask=None):
    """Find the cross edges of an image using the Roberts' Cross operator.

    The kernel is applied to the input image to produce separate measurements
    of the gradient component one orientation.

    Parameters
    ----------
    image : 2-D array
        Image to process.
    mask : 2-D array, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : 2-D array
        The Robert's edge map.

    Notes
    -----
    We use the following kernel::

      0   1
     -1   0

    """
    assert_nD(image, 2)
    image = img_as_float(image)
    result = convolve(image, ROBERTS_ND_WEIGHTS)
    return _mask_filter_result(result, mask)


def laplace(image, ksize=3, mask=None):
    """Find the edges of an image using the Laplace operator.

    Parameters
    ----------
    image : ndarray
        Image to process.
    ksize : int, optional
        Define the size of the discrete Laplacian operator such that it
        will have a size of (ksize,) * image.ndim.
    mask : ndarray, optional
        An optional mask to limit the application to a certain area.
        Note that pixels surrounding masked regions are also masked to
        prevent masked regions from affecting the result.

    Returns
    -------
    output : ndarray
        The Laplace edge map.

    Notes
    -----
    The Laplacian operator is generated using the function
    skimage.restoration.uft.laplacian().

    """
    image = img_as_float(image)
    # Create the discrete Laplacian operator - We keep only the real part of the filter
    _, laplace_op = laplacian(image.ndim, (ksize, ) * image.ndim)
    result = convolve(image, laplace_op)
    return _mask_filter_result(result, mask)
