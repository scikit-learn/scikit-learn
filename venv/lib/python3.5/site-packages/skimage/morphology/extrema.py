"""extrema.py - local minima and maxima

This module provides functions to find local maxima and minima of an image.
Here, local maxima (minima) are defined as connected sets of pixels with equal
grey level which is strictly greater (smaller) than the grey level of all
pixels in direct neighborhood of the connected set. In addition, the module
provides the closely related functions h-maxima and h-minima.

Soille, P. (2003). Morphological Image Analysis: Principles and Applications
(2nd ed.), Chapter 6. Springer-Verlag New York, Inc.
"""
from __future__ import absolute_import

import numpy as np

from . import greyreconstruct
from ..util import dtype_limits


def _add_constant_clip(image, const_value):
    """Add constant to the image while handling overflow issues gracefully.
    """
    min_dtype, max_dtype = dtype_limits(image, clip_negative=False)

    if const_value > (max_dtype - min_dtype):
        raise ValueError("The added constant is not compatible"
                         "with the image data type.")

    result = image + const_value
    result[image > max_dtype-const_value] = max_dtype
    return(result)


def _subtract_constant_clip(image, const_value):
    """Subtract constant from image while handling underflow issues.
    """
    min_dtype, max_dtype = dtype_limits(image, clip_negative=False)

    if const_value > (max_dtype-min_dtype):
        raise ValueError("The subtracted constant is not compatible"
                         "with the image data type.")

    result = image - const_value
    result[image < (const_value + min_dtype)] = min_dtype
    return(result)


def h_maxima(image, h, selem=None):
    """Determine all maxima of the image with height >= h.

    The local maxima are defined as connected sets of pixels with equal
    grey level strictly greater than the grey level of all pixels in direct
    neighborhood of the set.

    A local maximum M of height h is a local maximum for which
    there is at least one path joining M with a higher maximum on which the
    minimal value is f(M) - h (i.e. the values along the path are not
    decreasing by more than h with respect to the maximum's value) and no
    path for which the minimal value is greater.

    Parameters
    ----------
    image : ndarray
        The input image for which the maxima are to be calculated.
    h : unsigned integer
        The minimal height of all extracted maxima.
    selem : ndarray, optional
        The neighborhood expressed as an n-D array of 1's and 0's.
        Default is the ball of radius 1 according to the maximum norm
        (i.e. a 3x3 square for 2D images, a 3x3x3 cube for 3D images, etc.)

    Returns
    -------
    h_max : ndarray
       The maxima of height >= h. The result image is a binary image, where
       pixels belonging to the selected maxima take value 1, the others
       take value 0.

    See also
    --------
    skimage.morphology.extrema.h_minima
    skimage.morphology.extrema.local_maxima
    skimage.morphology.extrema.local_minima

    References
    ----------
    .. [1] Soille, P., "Morphological Image Analysis: Principles and
           Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.morphology import extrema

    We create an image (quadratic function with a maximum in the center and
    4 additional constant maxima.
    The heights of the maxima are: 1, 21, 41, 61, 81, 101

    >>> w = 10
    >>> x, y = np.mgrid[0:w,0:w]
    >>> f = 20 - 0.2*((x - w/2)**2 + (y-w/2)**2)
    >>> f[2:4,2:4] = 40; f[2:4,7:9] = 60; f[7:9,2:4] = 80; f[7:9,7:9] = 100
    >>> f = f.astype(np.int)

    We can calculate all maxima with a height of at least 40:

    >>> maxima = extrema.h_maxima(f, 40)

    The resulting image will contain 4 local maxima.
    """
    if np.issubdtype(image.dtype, np.floating):
        resolution = 2 * np.finfo(image.dtype).resolution
        if h < resolution:
            h = resolution
        h_corrected = h - resolution / 2.0
        shifted_img = image - h
    else:
        shifted_img = _subtract_constant_clip(image, h)
        h_corrected = h

    rec_img = greyreconstruct.reconstruction(shifted_img, image,
                                             method='dilation', selem=selem)
    residue_img = image - rec_img
    h_max = np.zeros(image.shape, dtype=np.uint8)
    h_max[residue_img >= h_corrected] = 1
    return h_max


def h_minima(image, h, selem=None):
    """Determine all minima of the image with depth >= h.

    The local minima are defined as connected sets of pixels with equal
    grey level strictly smaller than the grey levels of all pixels in direct
    neighborhood of the set.

    A local minimum M of depth h is a local minimum for which
    there is at least one path joining M with a deeper minimum on which the
    maximal value is f(M) + h (i.e. the values along the path are not
    increasing by more than h with respect to the minimum's value) and no
    path for which the maximal value is smaller.

    Parameters
    ----------
    image : ndarray
        The input image for which the minima are to be calculated.
    h : unsigned integer
        The minimal depth of all extracted minima.
    selem : ndarray, optional
        The neighborhood expressed as an n-D array of 1's and 0's.
        Default is the ball of radius 1 according to the maximum norm
        (i.e. a 3x3 square for 2D images, a 3x3x3 cube for 3D images, etc.)

    Returns
    -------
    h_min : ndarray
       The minima of depth >= h. The result image is a binary image, where
       pixels belonging to the selected minima take value 1, the other pixels
       take value 0.

    See also
    --------
    skimage.morphology.extrema.h_maxima
    skimage.morphology.extrema.local_maxima
    skimage.morphology.extrema.local_minima

    References
    ----------
    .. [1] Soille, P., "Morphological Image Analysis: Principles and
           Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.morphology import extrema

    We create an image (quadratic function with a minimum in the center and
    4 additional constant maxima.
    The depth of the minima are: 1, 21, 41, 61, 81, 101

    >>> w = 10
    >>> x, y = np.mgrid[0:w,0:w]
    >>> f = 180 + 0.2*((x - w/2)**2 + (y-w/2)**2)
    >>> f[2:4,2:4] = 160; f[2:4,7:9] = 140; f[7:9,2:4] = 120; f[7:9,7:9] = 100
    >>> f = f.astype(np.int)

    We can calculate all minima with a depth of at least 40:

    >>> minima = extrema.h_minima(f, 40)

    The resulting image will contain 4 local minima.
    """
    if np.issubdtype(image.dtype, np.floating):
        resolution = 2 * np.finfo(image.dtype).resolution
        if h < resolution:
            h = resolution
        h_corrected = h - resolution / 2.0
        shifted_img = image + h
    else:
        shifted_img = _add_constant_clip(image, h)
        h_corrected = h

    rec_img = greyreconstruct.reconstruction(shifted_img, image,
                                             method='erosion', selem=selem)
    residue_img = rec_img - image
    h_min = np.zeros(image.shape, dtype=np.uint8)
    h_min[residue_img >= h_corrected] = 1
    return h_min


def _find_min_diff(image):
    """
    Find the minimal difference of grey levels inside the image.
    """
    img_vec = np.unique(image.flatten())

    if img_vec.size == 1:
        return 0

    min_diff = np.min(img_vec[1:] - img_vec[:-1])
    return min_diff


def local_maxima(image, selem=None):
    """Determine all local maxima of the image.

    The local maxima are defined as connected sets of pixels with equal
    grey level strictly greater than the grey levels of all pixels in direct
    neighborhood of the set.

    For integer typed images, this corresponds to the h-maxima with h=1.
    For float typed images, h is determined as the smallest difference
    between grey levels.

    Parameters
    ----------
    image : ndarray
        The input image for which the maxima are to be calculated.
    selem : ndarray, optional
        The neighborhood expressed as an n-D array of 1's and 0's.
        Default is the ball of radius 1 according to the maximum norm
        (i.e. a 3x3 square for 2D images, a 3x3x3 cube for 3D images, etc.)

    Returns
    -------
    local_max : ndarray
       All local maxima of the image. The result image is a binary image,
       where pixels belonging to local maxima take value 1, the other pixels
       take value 0.

    See also
    --------
    skimage.morphology.extrema.h_minima
    skimage.morphology.extrema.h_maxima
    skimage.morphology.extrema.local_minima

    References
    ----------
    .. [1] Soille, P., "Morphological Image Analysis: Principles and
           Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.morphology import extrema

    We create an image (quadratic function with a maximum in the center and
    4 additional constant maxima.
    The heights of the maxima are: 1, 21, 41, 61, 81, 101

    >>> w = 10
    >>> x, y = np.mgrid[0:w,0:w]
    >>> f = 20 - 0.2*((x - w/2)**2 + (y-w/2)**2)
    >>> f[2:4,2:4] = 40; f[2:4,7:9] = 60; f[7:9,2:4] = 80; f[7:9,7:9] = 100
    >>> f = f.astype(np.int)

    We can calculate all local maxima:

    >>> maxima = extrema.local_maxima(f)

    The resulting image will contain all 6 local maxima.
    """
    # find the minimal grey level difference
    h = _find_min_diff(image)
    if h == 0:
        return np.zeros(image.shape, np.uint8)
    if not np.issubdtype(image.dtype, np.floating):
        h = 1
    local_max = h_maxima(image, h, selem=selem)
    return local_max


def local_minima(image, selem=None):
    """Determine all local minima of the image.

    The local minima are defined as connected sets of pixels with equal
    grey level strictly smaller than the grey levels of all pixels in direct
    neighborhood of the set.

    For integer typed images, this corresponds to the h-minima with h=1.
    For float typed images, h is determined as the smallest difference
    between grey levels.

    Parameters
    ----------
    image : ndarray
        The input image for which the minima are to be calculated.
    selem : ndarray, optional
        The neighborhood expressed as an n-D array of 1's and 0's.
        Default is the ball of radius 1 according to the maximum norm
        (i.e. a 3x3 square for 2D images, a 3x3x3 cube for 3D images, etc.)

    Returns
    -------
    local_min : ndarray
       All local minima of the image. The result image is a binary image,
       where pixels belonging to local minima take value 1, the other pixels
       take value 0.

    See also
    --------
    skimage.morphology.extrema.h_minima
    skimage.morphology.extrema.h_maxima
    skimage.morphology.extrema.local_maxima

    References
    ----------
    .. [1] Soille, P., "Morphological Image Analysis: Principles and
           Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.morphology import extrema

    We create an image (quadratic function with a minimum in the center and
    4 additional constant maxima.
    The depth of the minima are: 1, 21, 41, 61, 81, 101

    >>> w = 10
    >>> x, y = np.mgrid[0:w,0:w]
    >>> f = 180 + 0.2*((x - w/2)**2 + (y-w/2)**2)
    >>> f[2:4,2:4] = 160; f[2:4,7:9] = 140; f[7:9,2:4] = 120; f[7:9,7:9] = 100
    >>> f = f.astype(np.int)

    We can calculate all local minima:

    >>> minima = extrema.local_minima(f)

    The resulting image will contain all 6 local minima.
    """
    # find the minimal grey level difference
    h = _find_min_diff(image)
    if h == 0:
        return np.zeros(image.shape, np.uint8)
    if not np.issubdtype(image.dtype, np.floating):
        h = 1
    local_min = h_minima(image, h, selem=selem)
    return local_min
