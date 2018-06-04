#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Functions for converting between color spaces.

The "central" color space in this module is RGB, more specifically the linear
sRGB color space using D65 as a white-point [1]_.  This represents a
standard monitor (w/o gamma correction). For a good FAQ on color spaces see
[2]_.

The API consists of functions to convert to and from RGB as defined above, as
well as a generic function to convert to and from any supported color space
(which is done through RGB in most cases).


Supported color spaces
----------------------
* RGB : Red Green Blue.
        Here the sRGB standard [1]_.
* HSV : Hue, Saturation, Value.
        Uniquely defined when related to sRGB [3]_.
* RGB CIE : Red Green Blue.
        The original RGB CIE standard from 1931 [4]_. Primary colors are 700 nm
        (red), 546.1 nm (blue) and 435.8 nm (green).
* XYZ CIE : XYZ
        Derived from the RGB CIE color space. Chosen such that
        ``x == y == z == 1/3`` at the whitepoint, and all color matching
        functions are greater than zero everywhere.
* LAB CIE : Lightness, a, b
        Colorspace derived from XYZ CIE that is intended to be more
        perceptually uniform
* LUV CIE : Lightness, u, v
        Colorspace derived from XYZ CIE that is intended to be more
        perceptually uniform
* LCH CIE : Lightness, Chroma, Hue
        Defined in terms of LAB CIE.  C and H are the polar representation of
        a and b.  The polar angle C is defined to be on ``(0, 2*pi)``

:author: Nicolas Pinto (rgb2hsv)
:author: Ralf Gommers (hsv2rgb)
:author: Travis Oliphant (XYZ and RGB CIE functions)
:author: Matt Terry (lab2lch)
:author: Alex Izvorski (yuv2rgb, rgb2yuv and related)

:license: modified BSD

References
----------
.. [1] Official specification of sRGB, IEC 61966-2-1:1999.
.. [2] http://www.poynton.com/ColorFAQ.html
.. [3] http://en.wikipedia.org/wiki/HSL_and_HSV
.. [4] http://en.wikipedia.org/wiki/CIE_1931_color_space
"""

from __future__ import division

from warnings import warn
import numpy as np
from scipy import linalg
from ..util import dtype, dtype_limits


def guess_spatial_dimensions(image):
    """Make an educated guess about whether an image has a channels dimension.

    Parameters
    ----------
    image : ndarray
        The input image.

    Returns
    -------
    spatial_dims : int or None
        The number of spatial dimensions of `image`. If ambiguous, the value
        is ``None``.

    Raises
    ------
    ValueError
        If the image array has less than two or more than four dimensions.
    """
    if image.ndim == 2:
        return 2
    if image.ndim == 3 and image.shape[-1] != 3:
        return 3
    if image.ndim == 3 and image.shape[-1] == 3:
        return None
    if image.ndim == 4 and image.shape[-1] == 3:
        return 3
    else:
        raise ValueError("Expected 2D, 3D, or 4D array, got %iD." % image.ndim)


def convert_colorspace(arr, fromspace, tospace):
    """Convert an image array to a new color space.

    Valid color spaces are:
        'RGB', 'HSV', 'RGB CIE', 'XYZ', 'YUV', 'YIQ', 'YPbPr', 'YCbCr', 'YDbDr'

    Parameters
    ----------
    arr : array_like
        The image to convert.
    fromspace : valid color space
        The color space to convert from. Can be specified in lower case.
    tospace : valid color space
        The color space to convert to. Can be specified in lower case.

    Returns
    -------
    out : ndarray
        The converted image.

    Notes
    -----
    Conversion is performed through the "central" RGB color space,
    i.e. conversion from XYZ to HSV is implemented as ``XYZ -> RGB -> HSV``
    instead of directly.

    Examples
    --------
    >>> from skimage import data
    >>> img = data.astronaut()
    >>> img_hsv = convert_colorspace(img, 'RGB', 'HSV')

    """
    fromdict = {'rgb': lambda im: im, 'hsv': hsv2rgb, 'rgb cie': rgbcie2rgb,
                'xyz': xyz2rgb, 'yuv': yuv2rgb, 'yiq': yiq2rgb,
                'ypbpr': ypbpr2rgb, 'ycbcr': ycbcr2rgb, 'ydbdr': ydbdr2rgb}
    todict = {'rgb': lambda im: im, 'hsv': rgb2hsv, 'rgb cie': rgb2rgbcie,
              'xyz': rgb2xyz, 'yuv': rgb2yuv, 'yiq': rgb2yiq,
              'ypbpr': rgb2ypbpr, 'ycbcr': rgb2ycbcr, 'ydbdr': rgb2ydbdr}

    fromspace = fromspace.lower()
    tospace = tospace.lower()
    if fromspace not in fromdict:
        msg = '`fromspace` has to be one of {}'.format(fromdict.keys())
        raise ValueError(msg)
    if tospace not in todict:
        msg = '`tospace` has to be one of {}'.format(todict.keys())
        raise ValueError(msg)

    return todict[tospace](fromdict[fromspace](arr))


def _prepare_colorarray(arr):
    """Check the shape of the array and convert it to
    floating point representation.

    """
    arr = np.asanyarray(arr)

    if arr.ndim not in [3, 4] or arr.shape[-1] != 3:
        msg = ("the input array must be have a shape == (.., ..,[ ..,] 3)), " +
               "got (" + (", ".join(map(str, arr.shape))) + ")")
        raise ValueError(msg)

    return dtype.img_as_float(arr)


def _prepare_rgba_array(arr):
    """Check the shape of the array to be RGBA and convert it to
    floating point representation.

    """
    arr = np.asanyarray(arr)

    if arr.ndim not in [3, 4] or arr.shape[-1] != 4:
        msg = ("the input array must have a shape == (.., ..,[ ..,] 4)), "
               "got {0}".format(arr.shape))
        raise ValueError(msg)

    return dtype.img_as_float(arr)


def rgba2rgb(rgba, background=(1, 1, 1)):
    """RGBA to RGB conversion.

    Parameters
    ----------
    rgba : array_like
        The image in RGBA format, in a 3-D array of shape ``(.., .., 4)``.
    background : array_like
        The color of the background to blend the image with. A tuple
        containing 3 floats between 0 to 1 - the RGB value of the background.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `rgba` is not a 3-D array of shape ``(.., .., 4)``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Alpha_compositing#Alpha_blending

    Examples
    --------
    >>> from skimage import color
    >>> from skimage import data
    >>> img_rgba = data.logo()
    >>> img_rgb = color.rgba2rgb(img_rgba)
    """
    arr = _prepare_rgba_array(rgba)
    if isinstance(background, tuple) and len(background) != 3:
        raise ValueError('the background must be a tuple with 3 items - the '
                         'RGB color of the background. Got {0} items.'
                         .format(len(background)))

    alpha = arr[..., -1]
    channels = arr[..., :-1]
    out = np.empty_like(channels)

    for ichan in range(channels.shape[-1]):
        out[..., ichan] = np.clip(
            (1 - alpha) * background[ichan] + alpha * channels[..., ichan],
            a_min=0, a_max=1)
    return out


def rgb2hsv(rgb):
    """RGB to HSV color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Returns
    -------
    out : ndarray
        The image in HSV format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3-D array of shape ``(.., .., 3)``.

    Notes
    -----
    Conversion between RGB and HSV color spaces results in some loss of
    precision, due to integer arithmetic and rounding [1]_.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/HSL_and_HSV

    Examples
    --------
    >>> from skimage import color
    >>> from skimage import data
    >>> img = data.astronaut()
    >>> img_hsv = color.rgb2hsv(img)
    """
    arr = _prepare_colorarray(rgb)
    out = np.empty_like(arr)

    # -- V channel
    out_v = arr.max(-1)

    # -- S channel
    delta = arr.ptp(-1)
    # Ignore warning for zero divided by zero
    old_settings = np.seterr(invalid='ignore')
    out_s = delta / out_v
    out_s[delta == 0.] = 0.

    # -- H channel
    # red is max
    idx = (arr[:, :, 0] == out_v)
    out[idx, 0] = (arr[idx, 1] - arr[idx, 2]) / delta[idx]

    # green is max
    idx = (arr[:, :, 1] == out_v)
    out[idx, 0] = 2. + (arr[idx, 2] - arr[idx, 0]) / delta[idx]

    # blue is max
    idx = (arr[:, :, 2] == out_v)
    out[idx, 0] = 4. + (arr[idx, 0] - arr[idx, 1]) / delta[idx]
    out_h = (out[:, :, 0] / 6.) % 1.
    out_h[delta == 0.] = 0.

    np.seterr(**old_settings)

    # -- output
    out[:, :, 0] = out_h
    out[:, :, 1] = out_s
    out[:, :, 2] = out_v

    # remove NaN
    out[np.isnan(out)] = 0

    return out


def hsv2rgb(hsv):
    """HSV to RGB color space conversion.

    Parameters
    ----------
    hsv : array_like
        The image in HSV format, in a 3-D array of shape ``(.., .., 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `hsv` is not a 3-D array of shape ``(.., .., 3)``.

    Notes
    -----
    Conversion between RGB and HSV color spaces results in some loss of
    precision, due to integer arithmetic and rounding [1]_.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/HSL_and_HSV

    Examples
    --------
    >>> from skimage import data
    >>> img = data.astronaut()
    >>> img_hsv = rgb2hsv(img)
    >>> img_rgb = hsv2rgb(img_hsv)
    """
    arr = _prepare_colorarray(hsv)

    hi = np.floor(arr[:, :, 0] * 6)
    f = arr[:, :, 0] * 6 - hi
    p = arr[:, :, 2] * (1 - arr[:, :, 1])
    q = arr[:, :, 2] * (1 - f * arr[:, :, 1])
    t = arr[:, :, 2] * (1 - (1 - f) * arr[:, :, 1])
    v = arr[:, :, 2]

    hi = np.dstack([hi, hi, hi]).astype(np.uint8) % 6
    out = np.choose(hi, [np.dstack((v, t, p)),
                         np.dstack((q, v, p)),
                         np.dstack((p, v, t)),
                         np.dstack((p, q, v)),
                         np.dstack((t, p, v)),
                         np.dstack((v, p, q))])

    return out


# ---------------------------------------------------------------
# Primaries for the coordinate systems
# ---------------------------------------------------------------
cie_primaries = np.array([700, 546.1, 435.8])
sb_primaries = np.array([1. / 155, 1. / 190, 1. / 225]) * 1e5

# ---------------------------------------------------------------
# Matrices that define conversion between different color spaces
# ---------------------------------------------------------------

# From sRGB specification
xyz_from_rgb = np.array([[0.412453, 0.357580, 0.180423],
                         [0.212671, 0.715160, 0.072169],
                         [0.019334, 0.119193, 0.950227]])

rgb_from_xyz = linalg.inv(xyz_from_rgb)

# From http://en.wikipedia.org/wiki/CIE_1931_color_space
# Note: Travis's code did not have the divide by 0.17697
xyz_from_rgbcie = np.array([[0.49, 0.31, 0.20],
                            [0.17697, 0.81240, 0.01063],
                            [0.00, 0.01, 0.99]]) / 0.17697

rgbcie_from_xyz = linalg.inv(xyz_from_rgbcie)

# construct matrices to and from rgb:
rgbcie_from_rgb = np.dot(rgbcie_from_xyz, xyz_from_rgb)
rgb_from_rgbcie = np.dot(rgb_from_xyz, xyz_from_rgbcie)


gray_from_rgb = np.array([[0.2125, 0.7154, 0.0721],
                          [0, 0, 0],
                          [0, 0, 0]])

yuv_from_rgb = np.array([[ 0.299     ,  0.587     ,  0.114      ],
                         [-0.14714119, -0.28886916,  0.43601035 ],
                         [ 0.61497538, -0.51496512, -0.10001026 ]])

rgb_from_yuv = linalg.inv(yuv_from_rgb)

yiq_from_rgb = np.array([[0.299     ,  0.587     ,  0.114     ],
                         [0.59590059, -0.27455667, -0.32134392],
                         [0.21153661, -0.52273617,  0.31119955]])

rgb_from_yiq = linalg.inv(yiq_from_rgb)

ypbpr_from_rgb = np.array([[ 0.299   , 0.587   , 0.114   ],
                           [-0.168736,-0.331264, 0.5     ],
                           [ 0.5     ,-0.418688,-0.081312]])

rgb_from_ypbpr = linalg.inv(ypbpr_from_rgb)

ycbcr_from_rgb = np.array([[    65.481,   128.553,    24.966],
                           [   -37.797,   -74.203,   112.0  ],
                           [   112.0  ,   -93.786,   -18.214]])

rgb_from_ycbcr = linalg.inv(ycbcr_from_rgb)

ydbdr_from_rgb = np.array([[    0.299,   0.587,    0.114],
                           [   -0.45 ,  -0.883,    1.333],
                           [   -1.333,   1.116,    0.217]])

rgb_from_ydbdr = linalg.inv(ydbdr_from_rgb)


# CIE LAB constants for Observer=2A, Illuminant=D65
# NOTE: this is actually the XYZ values for the illuminant above.
lab_ref_white = np.array([0.95047, 1., 1.08883])

# XYZ coordinates of the illuminants, scaled to [0, 1]. For each illuminant I
# we have:
#
#   illuminant[I][0] corresponds to the XYZ coordinates for the 2 degree
#   field of view.
#
#   illuminant[I][1] corresponds to the XYZ coordinates for the 10 degree
#   field of view.
#
# The XYZ coordinates are calculated from [1], using the formula:
#
#   X = x * ( Y / y )
#   Y = Y
#   Z = ( 1 - x - y ) * ( Y / y )
#
# where Y = 1. The only exception is the illuminant "D65" with aperture angle
# 2, whose coordinates are copied from 'lab_ref_white' for
# backward-compatibility reasons.
#
#     References
#    ----------
#    .. [1] http://en.wikipedia.org/wiki/Standard_illuminant

illuminants = \
    {"A": {'2': (1.098466069456375, 1, 0.3558228003436005),
           '10': (1.111420406956693, 1, 0.3519978321919493)},
     "D50": {'2': (0.9642119944211994, 1, 0.8251882845188288),
             '10': (0.9672062750333777, 1, 0.8142801513128616)},
     "D55": {'2': (0.956797052643698, 1, 0.9214805860173273),
             '10': (0.9579665682254781, 1, 0.9092525159847462)},
     "D65": {'2': (0.95047, 1., 1.08883),   # This was: `lab_ref_white`
             '10': (0.94809667673716, 1, 1.0730513595166162)},
     "D75": {'2': (0.9497220898840717, 1, 1.226393520724154),
             '10': (0.9441713925645873, 1, 1.2064272211720228)},
     "E": {'2': (1.0, 1.0, 1.0),
           '10': (1.0, 1.0, 1.0)}}


def get_xyz_coords(illuminant, observer):
    """Get the XYZ coordinates of the given illuminant and observer [1]_.

    Parameters
    ----------
    illuminant : {"A", "D50", "D55", "D65", "D75", "E"}, optional
        The name of the illuminant (the function is NOT case sensitive).
    observer : {"2", "10"}, optional
        The aperture angle of the observer.

    Returns
    -------
    (x, y, z) : tuple
        A tuple with 3 elements containing the XYZ coordinates of the given
        illuminant.

    Raises
    ------
    ValueError
        If either the illuminant or the observer angle are not supported or
        unknown.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Standard_illuminant

    """
    illuminant = illuminant.upper()
    try:
        return illuminants[illuminant][observer]
    except KeyError:
        raise ValueError("Unknown illuminant/observer combination\
        (\'{0}\', \'{1}\')".format(illuminant, observer))

# Haematoxylin-Eosin-DAB colorspace
# From original Ruifrok's paper: A. C. Ruifrok and D. A. Johnston,
# "Quantification of histochemical staining by color deconvolution.,"
# Analytical and quantitative cytology and histology / the International
# Academy of Cytology [and] American Society of Cytology, vol. 23, no. 4,
# pp. 291-9, Aug. 2001.
rgb_from_hed = np.array([[0.65, 0.70, 0.29],
                         [0.07, 0.99, 0.11],
                         [0.27, 0.57, 0.78]])
hed_from_rgb = linalg.inv(rgb_from_hed)

# Following matrices are adapted form the Java code written by G.Landini.
# The original code is available at:
# http://www.dentistry.bham.ac.uk/landinig/software/cdeconv/cdeconv.html

# Hematoxylin + DAB
rgb_from_hdx = np.array([[0.650, 0.704, 0.286],
                         [0.268, 0.570, 0.776],
                         [0.0, 0.0, 0.0]])
rgb_from_hdx[2, :] = np.cross(rgb_from_hdx[0, :], rgb_from_hdx[1, :])
hdx_from_rgb = linalg.inv(rgb_from_hdx)

# Feulgen + Light Green
rgb_from_fgx = np.array([[0.46420921, 0.83008335, 0.30827187],
                         [0.94705542, 0.25373821, 0.19650764],
                         [0.0, 0.0, 0.0]])
rgb_from_fgx[2, :] = np.cross(rgb_from_fgx[0, :], rgb_from_fgx[1, :])
fgx_from_rgb = linalg.inv(rgb_from_fgx)

# Giemsa: Methyl Blue + Eosin
rgb_from_bex = np.array([[0.834750233, 0.513556283, 0.196330403],
                         [0.092789, 0.954111, 0.283111],
                         [0.0, 0.0, 0.0]])
rgb_from_bex[2, :] = np.cross(rgb_from_bex[0, :], rgb_from_bex[1, :])
bex_from_rgb = linalg.inv(rgb_from_bex)

# FastRed + FastBlue +  DAB
rgb_from_rbd = np.array([[0.21393921, 0.85112669, 0.47794022],
                         [0.74890292, 0.60624161, 0.26731082],
                         [0.268, 0.570, 0.776]])
rbd_from_rgb = linalg.inv(rgb_from_rbd)

# Methyl Green + DAB
rgb_from_gdx = np.array([[0.98003, 0.144316, 0.133146],
                         [0.268, 0.570, 0.776],
                         [0.0, 0.0, 0.0]])
rgb_from_gdx[2, :] = np.cross(rgb_from_gdx[0, :], rgb_from_gdx[1, :])
gdx_from_rgb = linalg.inv(rgb_from_gdx)

# Hematoxylin + AEC
rgb_from_hax = np.array([[0.650, 0.704, 0.286],
                         [0.2743, 0.6796, 0.6803],
                         [0.0, 0.0, 0.0]])
rgb_from_hax[2, :] = np.cross(rgb_from_hax[0, :], rgb_from_hax[1, :])
hax_from_rgb = linalg.inv(rgb_from_hax)

# Blue matrix Anilline Blue + Red matrix Azocarmine + Orange matrix Orange-G
rgb_from_bro = np.array([[0.853033, 0.508733, 0.112656],
                         [0.09289875, 0.8662008, 0.49098468],
                         [0.10732849, 0.36765403, 0.9237484]])
bro_from_rgb = linalg.inv(rgb_from_bro)

# Methyl Blue + Ponceau Fuchsin
rgb_from_bpx = np.array([[0.7995107, 0.5913521, 0.10528667],
                         [0.09997159, 0.73738605, 0.6680326],
                         [0.0, 0.0, 0.0]])
rgb_from_bpx[2, :] = np.cross(rgb_from_bpx[0, :], rgb_from_bpx[1, :])
bpx_from_rgb = linalg.inv(rgb_from_bpx)

# Alcian Blue + Hematoxylin
rgb_from_ahx = np.array([[0.874622, 0.457711, 0.158256],
                         [0.552556, 0.7544, 0.353744],
                         [0.0, 0.0, 0.0]])
rgb_from_ahx[2, :] = np.cross(rgb_from_ahx[0, :], rgb_from_ahx[1, :])
ahx_from_rgb = linalg.inv(rgb_from_ahx)

# Hematoxylin + PAS
rgb_from_hpx = np.array([[0.644211, 0.716556, 0.266844],
                         [0.175411, 0.972178, 0.154589],
                         [0.0, 0.0, 0.0]])
rgb_from_hpx[2, :] = np.cross(rgb_from_hpx[0, :], rgb_from_hpx[1, :])
hpx_from_rgb = linalg.inv(rgb_from_hpx)

# -------------------------------------------------------------
# The conversion functions that make use of the matrices above
# -------------------------------------------------------------


def _convert(matrix, arr):
    """Do the color space conversion.

    Parameters
    ----------
    matrix : array_like
        The 3x3 matrix to use.
    arr : array_like
        The input array.

    Returns
    -------
    out : ndarray, dtype=float
        The converted array.
    """
    arr = _prepare_colorarray(arr)

    return np.dot(arr, matrix.T.copy())


def xyz2rgb(xyz):
    """XYZ to RGB color space conversion.

    Parameters
    ----------
    xyz : array_like
        The image in XYZ format, in a 3-D array of shape ``(.., .., 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `xyz` is not a 3-D array of shape ``(.., .., 3)``.

    Notes
    -----
    The CIE XYZ color space is derived from the CIE RGB color space. Note
    however that this function converts to sRGB.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/CIE_1931_color_space

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2xyz, xyz2rgb
    >>> img = data.astronaut()
    >>> img_xyz = rgb2xyz(img)
    >>> img_rgb = xyz2rgb(img_xyz)
    """
    # Follow the algorithm from http://www.easyrgb.com/index.php
    # except we don't multiply/divide by 100 in the conversion
    arr = _convert(rgb_from_xyz, xyz)
    mask = arr > 0.0031308
    arr[mask] = 1.055 * np.power(arr[mask], 1 / 2.4) - 0.055
    arr[~mask] *= 12.92
    arr[arr < 0] = 0
    arr[arr > 1] = 1
    return arr


def rgb2xyz(rgb):
    """RGB to XYZ color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3- or 4-D array of shape
        ``(.., ..,[ ..,] 3)``.

    Returns
    -------
    out : ndarray
        The image in XYZ format, in a 3- or 4-D array of shape
        ``(.., ..,[ ..,] 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3- or 4-D array of shape ``(.., ..,[ ..,] 3)``.

    Notes
    -----
    The CIE XYZ color space is derived from the CIE RGB color space. Note
    however that this function converts from sRGB.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/CIE_1931_color_space

    Examples
    --------
    >>> from skimage import data
    >>> img = data.astronaut()
    >>> img_xyz = rgb2xyz(img)
    """
    # Follow the algorithm from http://www.easyrgb.com/index.php
    # except we don't multiply/divide by 100 in the conversion
    arr = _prepare_colorarray(rgb).copy()
    mask = arr > 0.04045
    arr[mask] = np.power((arr[mask] + 0.055) / 1.055, 2.4)
    arr[~mask] /= 12.92
    return _convert(xyz_from_rgb, arr)


def rgb2rgbcie(rgb):
    """RGB to RGB CIE color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB CIE format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3-D array of shape ``(.., .., 3)``.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/CIE_1931_color_space

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2rgbcie
    >>> img = data.astronaut()
    >>> img_rgbcie = rgb2rgbcie(img)
    """
    return _convert(rgbcie_from_rgb, rgb)


def rgbcie2rgb(rgbcie):
    """RGB CIE to RGB color space conversion.

    Parameters
    ----------
    rgbcie : array_like
        The image in RGB CIE format, in a 3-D array of shape ``(.., .., 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `rgbcie` is not a 3-D array of shape ``(.., .., 3)``.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/CIE_1931_color_space

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2rgbcie, rgbcie2rgb
    >>> img = data.astronaut()
    >>> img_rgbcie = rgb2rgbcie(img)
    >>> img_rgb = rgbcie2rgb(img_rgbcie)
    """
    return _convert(rgb_from_rgbcie, rgbcie)


def rgb2gray(rgb):
    """Compute luminance of an RGB image.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3-D or 4-D array of shape
        ``(.., ..,[ ..,] 3)``, or in RGBA format with shape
        ``(.., ..,[ ..,] 4)``.

    Returns
    -------
    out : ndarray
        The luminance image - an array which is the same size as the input
        array, but with the channel dimension removed.

    Raises
    ------
    ValueError
        If `rgb2gray` is not a 3-D or 4-D arrays of shape
        ``(.., ..,[ ..,] 3)`` or ``(.., ..,[ ..,] 4)``.

    References
    ----------
    .. [1] http://www.poynton.com/PDFs/ColorFAQ.pdf

    Notes
    -----
    The weights used in this conversion are calibrated for contemporary
    CRT phosphors::

        Y = 0.2125 R + 0.7154 G + 0.0721 B

    If there is an alpha channel present, it is ignored.

    Examples
    --------
    >>> from skimage.color import rgb2gray
    >>> from skimage import data
    >>> img = data.astronaut()
    >>> img_gray = rgb2gray(img)
    """

    if rgb.ndim == 2:
        return np.ascontiguousarray(rgb)

    rgb = _prepare_colorarray(rgb[..., :3])

    gray = 0.2125 * rgb[..., 0]
    gray[:] += 0.7154 * rgb[..., 1]
    gray[:] += 0.0721 * rgb[..., 2]

    return gray


rgb2grey = rgb2gray


def gray2rgb(image, alpha=None):
    """Create an RGB representation of a gray-level image.

    Parameters
    ----------
    image : array_like
        Input image of shape ``(M[, N][, P])``.
    alpha : bool, optional
        Ensure that the output image has an alpha layer.  If None,
        alpha layers are passed through but not created.

    Returns
    -------
    rgb : ndarray
        RGB image of shape ``(M[, N][, P], 3)``.

    Raises
    ------
    ValueError
        If the input is not a 1-, 2- or 3-dimensional image.

    Notes
    -----
    If the input is a 1-dimensional image of shape ``(M, )``, the output
    will be shape ``(M, 3)``.
    """
    is_rgb = False
    is_alpha = False
    dims = np.squeeze(image).ndim

    if dims == 3:
        if image.shape[2] == 3:
            is_rgb = True
        elif image.shape[2] == 4:
            is_alpha = True
            is_rgb = True

    if is_rgb:
        if alpha is False:
            image = image[..., :3]

        elif alpha is True and is_alpha is False:
            alpha_layer = (np.ones_like(image[..., 0, np.newaxis]) *
                           dtype_limits(image, clip_negative=False)[1])
            image = np.concatenate((image, alpha_layer), axis=2)

        return image

    elif dims in (1, 2, 3):
        image = image[..., np.newaxis]

        if alpha:
            alpha_layer = (np.ones_like(image) * dtype_limits(image, clip_negative=False)[1])
            return np.concatenate(3 * (image,) + (alpha_layer,), axis=-1)
        else:
            return np.concatenate(3 * (image,), axis=-1)

    else:
        raise ValueError("Input image expected to be RGB, RGBA or gray.")

grey2rgb = gray2rgb


def xyz2lab(xyz, illuminant="D65", observer="2"):
    """XYZ to CIE-LAB color space conversion.

    Parameters
    ----------
    xyz : array_like
        The image in XYZ format, in a 3- or 4-D array of shape
        ``(.., ..,[ ..,] 3)``.
    illuminant : {"A", "D50", "D55", "D65", "D75", "E"}, optional
        The name of the illuminant (the function is NOT case sensitive).
    observer : {"2", "10"}, optional
        The aperture angle of the observer.

    Returns
    -------
    out : ndarray
        The image in CIE-LAB format, in a 3- or 4-D array of shape
        ``(.., ..,[ ..,] 3)``.

    Raises
    ------
    ValueError
        If `xyz` is not a 3-D array of shape ``(.., ..,[ ..,] 3)``.
    ValueError
        If either the illuminant or the observer angle is unsupported or
        unknown.

    Notes
    -----
    By default Observer= 2A, Illuminant= D65. CIE XYZ tristimulus values
    x_ref=95.047, y_ref=100., z_ref=108.883. See function `get_xyz_coords` for
    a list of supported illuminants.

    References
    ----------
    .. [1] http://www.easyrgb.com/index.php?X=MATH&H=07#text7
    .. [2] http://en.wikipedia.org/wiki/Lab_color_space

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2xyz, xyz2lab
    >>> img = data.astronaut()
    >>> img_xyz = rgb2xyz(img)
    >>> img_lab = xyz2lab(img_xyz)
    """
    arr = _prepare_colorarray(xyz)

    xyz_ref_white = get_xyz_coords(illuminant, observer)

    # scale by CIE XYZ tristimulus values of the reference white point
    arr = arr / xyz_ref_white

    # Nonlinear distortion and linear transformation
    mask = arr > 0.008856
    arr[mask] = np.power(arr[mask], 1. / 3.)
    arr[~mask] = 7.787 * arr[~mask] + 16. / 116.

    x, y, z = arr[..., 0], arr[..., 1], arr[..., 2]

    # Vector scaling
    L = (116. * y) - 16.
    a = 500.0 * (x - y)
    b = 200.0 * (y - z)

    return np.concatenate([x[..., np.newaxis] for x in [L, a, b]], axis=-1)


def lab2xyz(lab, illuminant="D65", observer="2"):
    """CIE-LAB to XYZcolor space conversion.

    Parameters
    ----------
    lab : array_like
        The image in lab format, in a 3-D array of shape ``(.., .., 3)``.
    illuminant : {"A", "D50", "D55", "D65", "D75", "E"}, optional
        The name of the illuminant (the function is NOT case sensitive).
    observer : {"2", "10"}, optional
        The aperture angle of the observer.

    Returns
    -------
    out : ndarray
        The image in XYZ format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `lab` is not a 3-D array of shape ``(.., .., 3)``.
    ValueError
        If either the illuminant or the observer angle are not supported or
        unknown.
    UserWarning
        If any of the pixels are invalid (Z < 0).


    Notes
    -----
    By default Observer= 2A, Illuminant= D65. CIE XYZ tristimulus values x_ref
    = 95.047, y_ref = 100., z_ref = 108.883. See function 'get_xyz_coords' for
    a list of supported illuminants.

    References
    ----------
    .. [1] http://www.easyrgb.com/index.php?X=MATH&H=07#text7
    .. [2] http://en.wikipedia.org/wiki/Lab_color_space

    """

    arr = _prepare_colorarray(lab).copy()

    L, a, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    y = (L + 16.) / 116.
    x = (a / 500.) + y
    z = y - (b / 200.)

    if np.any(z < 0):
        invalid = np.nonzero(z < 0)
        warn('Color data out of range: Z < 0 in %s pixels' % invalid[0].size)
        z[invalid] = 0

    out = np.dstack([x, y, z])

    mask = out > 0.2068966
    out[mask] = np.power(out[mask], 3.)
    out[~mask] = (out[~mask] - 16.0 / 116.) / 7.787

    # rescale to the reference white (illuminant)
    xyz_ref_white = get_xyz_coords(illuminant, observer)
    out *= xyz_ref_white
    return out


def rgb2lab(rgb, illuminant="D65", observer="2"):
    """RGB to lab color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3- or 4-D array of shape
        ``(.., ..,[ ..,] 3)``.
    illuminant : {"A", "D50", "D55", "D65", "D75", "E"}, optional
        The name of the illuminant (the function is NOT case sensitive).
    observer : {"2", "10"}, optional
        The aperture angle of the observer.

    Returns
    -------
    out : ndarray
        The image in Lab format, in a 3- or 4-D array of shape
        ``(.., ..,[ ..,] 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3- or 4-D array of shape ``(.., ..,[ ..,] 3)``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Standard_illuminant

    Notes
    -----
    This function uses rgb2xyz and xyz2lab.
    By default Observer= 2A, Illuminant= D65. CIE XYZ tristimulus values
    x_ref=95.047, y_ref=100., z_ref=108.883. See function `get_xyz_coords` for
    a list of supported illuminants.
    """
    return xyz2lab(rgb2xyz(rgb), illuminant, observer)


def lab2rgb(lab, illuminant="D65", observer="2"):
    """Lab to RGB color space conversion.

    Parameters
    ----------
    lab : array_like
        The image in Lab format, in a 3-D array of shape ``(.., .., 3)``.
    illuminant : {"A", "D50", "D55", "D65", "D75", "E"}, optional
        The name of the illuminant (the function is NOT case sensitive).
    observer : {"2", "10"}, optional
        The aperture angle of the observer.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `lab` is not a 3-D array of shape ``(.., .., 3)``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Standard_illuminant

    Notes
    -----
    This function uses lab2xyz and xyz2rgb.
    By default Observer= 2A, Illuminant= D65. CIE XYZ tristimulus values
    x_ref=95.047, y_ref=100., z_ref=108.883. See function `get_xyz_coords` for
    a list of supported illuminants.
    """
    return xyz2rgb(lab2xyz(lab, illuminant, observer))


def xyz2luv(xyz, illuminant="D65", observer="2"):
    """XYZ to CIE-Luv color space conversion.

    Parameters
    ----------
    xyz : (M, N, [P,] 3) array_like
        The 3 or 4 dimensional image in XYZ format. Final dimension denotes
        channels.
    illuminant : {"A", "D50", "D55", "D65", "D75", "E"}, optional
        The name of the illuminant (the function is NOT case sensitive).
    observer : {"2", "10"}, optional
        The aperture angle of the observer.

    Returns
    -------
    out : (M, N, [P,] 3) ndarray
        The image in CIE-Luv format. Same dimensions as input.

    Raises
    ------
    ValueError
        If `xyz` is not a 3-D or 4-D array of shape ``(M, N, [P,] 3)``.
    ValueError
        If either the illuminant or the observer angle are not supported or
        unknown.

    Notes
    -----
    By default XYZ conversion weights use observer=2A. Reference whitepoint
    for D65 Illuminant, with XYZ tristimulus values of ``(95.047, 100.,
    108.883)``. See function 'get_xyz_coords' for a list of supported
    illuminants.

    References
    ----------
    .. [1] http://www.easyrgb.com/index.php?X=MATH&H=16#text16
    .. [2] http://en.wikipedia.org/wiki/CIELUV

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2xyz, xyz2luv
    >>> img = data.astronaut()
    >>> img_xyz = rgb2xyz(img)
    >>> img_luv = xyz2luv(img_xyz)
    """
    arr = _prepare_colorarray(xyz)

    # extract channels
    x, y, z = arr[..., 0], arr[..., 1], arr[..., 2]

    eps = np.finfo(np.float).eps

    # compute y_r and L
    xyz_ref_white = get_xyz_coords(illuminant, observer)
    L = y / xyz_ref_white[1]
    mask = L > 0.008856
    L[mask] = 116. * np.power(L[mask], 1. / 3.) - 16.
    L[~mask] = 903.3 * L[~mask]

    u0 = 4 * xyz_ref_white[0] / np.dot([1, 15, 3], xyz_ref_white)
    v0 = 9 * xyz_ref_white[1] / np.dot([1, 15, 3], xyz_ref_white)

    # u' and v' helper functions
    def fu(X, Y, Z):
        return (4. * X) / (X + 15. * Y + 3. * Z + eps)

    def fv(X, Y, Z):
        return (9. * Y) / (X + 15. * Y + 3. * Z + eps)

    # compute u and v using helper functions
    u = 13. * L * (fu(x, y, z) - u0)
    v = 13. * L * (fv(x, y, z) - v0)

    return np.concatenate([q[..., np.newaxis] for q in [L, u, v]], axis=-1)


def luv2xyz(luv, illuminant="D65", observer="2"):
    """CIE-Luv to XYZ color space conversion.

    Parameters
    ----------
    luv : (M, N, [P,] 3) array_like
        The 3 or 4 dimensional image in CIE-Luv format. Final dimension denotes
        channels.
    illuminant : {"A", "D50", "D55", "D65", "D75", "E"}, optional
        The name of the illuminant (the function is NOT case sensitive).
    observer : {"2", "10"}, optional
        The aperture angle of the observer.

    Returns
    -------
    out : (M, N, [P,] 3) ndarray
        The image in XYZ format. Same dimensions as input.

    Raises
    ------
    ValueError
        If `luv` is not a 3-D or 4-D array of shape ``(M, N, [P,] 3)``.
    ValueError
        If either the illuminant or the observer angle are not supported or
        unknown.

    Notes
    -----
    XYZ conversion weights use observer=2A. Reference whitepoint for D65
    Illuminant, with XYZ tristimulus values of ``(95.047, 100., 108.883)``. See
    function 'get_xyz_coords' for a list of supported illuminants.

    References
    ----------
    .. [1] http://www.easyrgb.com/index.php?X=MATH&H=16#text16
    .. [2] http://en.wikipedia.org/wiki/CIELUV

    """

    arr = _prepare_colorarray(luv).copy()

    L, u, v = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]

    eps = np.finfo(np.float).eps

    # compute y
    y = L.copy()
    mask = y > 7.999625
    y[mask] = np.power((y[mask] + 16.) / 116., 3.)
    y[~mask] = y[~mask] / 903.3
    xyz_ref_white = get_xyz_coords(illuminant, observer)
    y *= xyz_ref_white[1]

    # reference white x,z
    uv_weights = [1, 15, 3]
    u0 = 4 * xyz_ref_white[0] / np.dot(uv_weights, xyz_ref_white)
    v0 = 9 * xyz_ref_white[1] / np.dot(uv_weights, xyz_ref_white)

    # compute intermediate values
    a = u0 + u / (13. * L + eps)
    b = v0 + v / (13. * L + eps)
    c = 3 * y * (5 * b - 3)

    # compute x and z
    z = ((a - 4) * c - 15 * a * b * y) / (12 * b)
    x = -(c / b + 3. * z)

    return np.concatenate([q[..., np.newaxis] for q in [x, y, z]], axis=-1)


def rgb2luv(rgb):
    """RGB to CIE-Luv color space conversion.

    Parameters
    ----------
    rgb : (M, N, [P,] 3) array_like
        The 3 or 4 dimensional image in RGB format. Final dimension denotes
        channels.

    Returns
    -------
    out : (M, N, [P,] 3) ndarray
        The image in CIE Luv format. Same dimensions as input.

    Raises
    ------
    ValueError
        If `rgb` is not a 3-D or 4-D array of shape ``(M, N, [P,] 3)``.

    Notes
    -----
    This function uses rgb2xyz and xyz2luv.

    References
    ----------
    .. [1] http://www.easyrgb.com/index.php?X=MATH&H=16#text16
    .. [2] http://www.easyrgb.com/index.php?X=MATH&H=02#text2
    .. [3] http://en.wikipedia.org/wiki/CIELUV

    """
    return xyz2luv(rgb2xyz(rgb))


def luv2rgb(luv):
    """Luv to RGB color space conversion.

    Parameters
    ----------
    luv : (M, N, [P,] 3) array_like
        The 3 or 4 dimensional image in CIE Luv format. Final dimension denotes
        channels.

    Returns
    -------
    out : (M, N, [P,] 3) ndarray
        The image in RGB format. Same dimensions as input.

    Raises
    ------
    ValueError
        If `luv` is not a 3-D or 4-D array of shape ``(M, N, [P,] 3)``.

    Notes
    -----
    This function uses luv2xyz and xyz2rgb.
    """
    return xyz2rgb(luv2xyz(luv))


def rgb2hed(rgb):
    """RGB to Haematoxylin-Eosin-DAB (HED) color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Returns
    -------
    out : ndarray
        The image in HED format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3-D array of shape ``(.., .., 3)``.


    References
    ----------
    .. [1] A. C. Ruifrok and D. A. Johnston, "Quantification of histochemical
           staining by color deconvolution.," Analytical and quantitative
           cytology and histology / the International Academy of Cytology [and]
           American Society of Cytology, vol. 23, no. 4, pp. 291-9, Aug. 2001.

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2hed
    >>> ihc = data.immunohistochemistry()
    >>> ihc_hed = rgb2hed(ihc)
    """
    return separate_stains(rgb, hed_from_rgb)


def hed2rgb(hed):
    """Haematoxylin-Eosin-DAB (HED) to RGB color space conversion.

    Parameters
    ----------
    hed : array_like
        The image in the HED color space, in a 3-D array of shape
        ``(.., .., 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `hed` is not a 3-D array of shape ``(.., .., 3)``.

    References
    ----------
    .. [1] A. C. Ruifrok and D. A. Johnston, "Quantification of histochemical
           staining by color deconvolution.," Analytical and quantitative
           cytology and histology / the International Academy of Cytology [and]
           American Society of Cytology, vol. 23, no. 4, pp. 291-9, Aug. 2001.

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2hed, hed2rgb
    >>> ihc = data.immunohistochemistry()
    >>> ihc_hed = rgb2hed(ihc)
    >>> ihc_rgb = hed2rgb(ihc_hed)
    """
    return combine_stains(hed, rgb_from_hed)


def separate_stains(rgb, conv_matrix):
    """RGB to stain color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.
    conv_matrix: ndarray
        The stain separation matrix as described by G. Landini [1]_.

    Returns
    -------
    out : ndarray
        The image in stain color space, in a 3-D array of shape
        ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3-D array of shape ``(.., .., 3)``.

    Notes
    -----
    Stain separation matrices available in the ``color`` module and their
    respective colorspace:

    * ``hed_from_rgb``: Hematoxylin + Eosin + DAB
    * ``hdx_from_rgb``: Hematoxylin + DAB
    * ``fgx_from_rgb``: Feulgen + Light Green
    * ``bex_from_rgb``: Giemsa stain : Methyl Blue + Eosin
    * ``rbd_from_rgb``: FastRed + FastBlue +  DAB
    * ``gdx_from_rgb``: Methyl Green + DAB
    * ``hax_from_rgb``: Hematoxylin + AEC
    * ``bro_from_rgb``: Blue matrix Anilline Blue + Red matrix Azocarmine\
                        + Orange matrix Orange-G
    * ``bpx_from_rgb``: Methyl Blue + Ponceau Fuchsin
    * ``ahx_from_rgb``: Alcian Blue + Hematoxylin
    * ``hpx_from_rgb``: Hematoxylin + PAS

    References
    ----------
    .. [1] http://www.dentistry.bham.ac.uk/landinig/software/cdeconv/cdeconv.html

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import separate_stains, hdx_from_rgb
    >>> ihc = data.immunohistochemistry()
    >>> ihc_hdx = separate_stains(ihc, hdx_from_rgb)
    """
    rgb = dtype.img_as_float(rgb, force_copy=True)
    rgb += 2
    stains = np.dot(np.reshape(-np.log(rgb), (-1, 3)), conv_matrix)
    return np.reshape(stains, rgb.shape)


def combine_stains(stains, conv_matrix):
    """Stain to RGB color space conversion.

    Parameters
    ----------
    stains : array_like
        The image in stain color space, in a 3-D array of shape
        ``(.., .., 3)``.
    conv_matrix: ndarray
        The stain separation matrix as described by G. Landini [1]_.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3-D array of shape ``(.., .., 3)``.

    Raises
    ------
    ValueError
        If `stains` is not a 3-D array of shape ``(.., .., 3)``.

    Notes
    -----
    Stain combination matrices available in the ``color`` module and their
    respective colorspace:

    * ``rgb_from_hed``: Hematoxylin + Eosin + DAB
    * ``rgb_from_hdx``: Hematoxylin + DAB
    * ``rgb_from_fgx``: Feulgen + Light Green
    * ``rgb_from_bex``: Giemsa stain : Methyl Blue + Eosin
    * ``rgb_from_rbd``: FastRed + FastBlue +  DAB
    * ``rgb_from_gdx``: Methyl Green + DAB
    * ``rgb_from_hax``: Hematoxylin + AEC
    * ``rgb_from_bro``: Blue matrix Anilline Blue + Red matrix Azocarmine\
                        + Orange matrix Orange-G
    * ``rgb_from_bpx``: Methyl Blue + Ponceau Fuchsin
    * ``rgb_from_ahx``: Alcian Blue + Hematoxylin
    * ``rgb_from_hpx``: Hematoxylin + PAS

    References
    ----------
    .. [1] http://www.dentistry.bham.ac.uk/landinig/software/cdeconv/cdeconv.html


    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import (separate_stains, combine_stains,
    ...                            hdx_from_rgb, rgb_from_hdx)
    >>> ihc = data.immunohistochemistry()
    >>> ihc_hdx = separate_stains(ihc, hdx_from_rgb)
    >>> ihc_rgb = combine_stains(ihc_hdx, rgb_from_hdx)
    """
    from ..exposure import rescale_intensity

    stains = dtype.img_as_float(stains)
    logrgb2 = np.dot(-np.reshape(stains, (-1, 3)), conv_matrix)
    rgb2 = np.exp(logrgb2)
    return rescale_intensity(np.reshape(rgb2 - 2, stains.shape),
                             in_range=(-1, 1))


def lab2lch(lab):
    """CIE-LAB to CIE-LCH color space conversion.

    LCH is the cylindrical representation of the LAB (Cartesian) colorspace

    Parameters
    ----------
    lab : array_like
        The N-D image in CIE-LAB format. The last (``N+1``-th) dimension must
        have at least 3 elements, corresponding to the ``L``, ``a``, and ``b``
        color channels.  Subsequent elements are copied.

    Returns
    -------
    out : ndarray
        The image in LCH format, in a N-D array with same shape as input `lab`.

    Raises
    ------
    ValueError
        If `lch` does not have at least 3 color channels (i.e. l, a, b).

    Notes
    -----
    The Hue is expressed as an angle between ``(0, 2*pi)``

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2lab, lab2lch
    >>> img = data.astronaut()
    >>> img_lab = rgb2lab(img)
    >>> img_lch = lab2lch(img_lab)
    """
    lch = _prepare_lab_array(lab)

    a, b = lch[..., 1], lch[..., 2]
    lch[..., 1], lch[..., 2] = _cart2polar_2pi(a, b)
    return lch


def _cart2polar_2pi(x, y):
    """convert cartesian coordinates to polar (uses non-standard theta range!)

    NON-STANDARD RANGE! Maps to ``(0, 2*pi)`` rather than usual ``(-pi, +pi)``
    """
    r, t = np.hypot(x, y), np.arctan2(y, x)
    t += np.where(t < 0., 2 * np.pi, 0)
    return r, t


def lch2lab(lch):
    """CIE-LCH to CIE-LAB color space conversion.

    LCH is the cylindrical representation of the LAB (Cartesian) colorspace

    Parameters
    ----------
    lch : array_like
        The N-D image in CIE-LCH format. The last (``N+1``-th) dimension must
        have at least 3 elements, corresponding to the ``L``, ``a``, and ``b``
        color channels.  Subsequent elements are copied.

    Returns
    -------
    out : ndarray
        The image in LAB format, with same shape as input `lch`.

    Raises
    ------
    ValueError
        If `lch` does not have at least 3 color channels (i.e. l, c, h).

    Examples
    --------
    >>> from skimage import data
    >>> from skimage.color import rgb2lab, lch2lab
    >>> img = data.astronaut()
    >>> img_lab = rgb2lab(img)
    >>> img_lch = lab2lch(img_lab)
    >>> img_lab2 = lch2lab(img_lch)
    """
    lch = _prepare_lab_array(lch)

    c, h = lch[..., 1], lch[..., 2]
    lch[..., 1], lch[..., 2] = c * np.cos(h), c * np.sin(h)
    return lch


def _prepare_lab_array(arr):
    """Ensure input for lab2lch, lch2lab are well-posed.

    Arrays must be in floating point and have at least 3 elements in
    last dimension.  Return a new array.
    """
    arr = np.asarray(arr)
    shape = arr.shape
    if shape[-1] < 3:
        raise ValueError('Input array has less than 3 color channels')
    return dtype.img_as_float(arr, force_copy=True)


def rgb2yuv(rgb):
    """RGB to YUV color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in YUV format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    Notes
    -----
    Y is between 0 and 1.  Use YCbCr instead of YUV for the color space which
    is commonly used by video codecs (where Y ranges from 16 to 235)

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YUV

    """
    return _convert(yuv_from_rgb, rgb)


def rgb2yiq(rgb):
    """RGB to YIQ color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in YIQ format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.
    """
    return _convert(yiq_from_rgb, rgb)


def rgb2ypbpr(rgb):
    """RGB to YPbPr color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in YPbPr format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YPbPr

    """
    return _convert(ypbpr_from_rgb, rgb)


def rgb2ycbcr(rgb):
    """RGB to YCbCr color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in YCbCr format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    Notes
    -----
    Y is between 16 and 235.  This is the color space which is commonly used
    by video codecs, it is sometimes incorrectly called "YUV"

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YCbCr

    """
    arr = _convert(ycbcr_from_rgb, rgb)
    arr[..., 0] += 16
    arr[..., 1] += 128
    arr[..., 2] += 128
    return arr


def rgb2ydbdr(rgb):
    """RGB to YDbDr color space conversion.

    Parameters
    ----------
    rgb : array_like
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in YDbDr format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `rgb` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    Notes
    -----
    This is the color space which is commonly used
    by video codecs, it is also the reversible color transform in JPEG2000.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YDbDr

    """
    arr = _convert(ydbdr_from_rgb, rgb)
    return arr


def yuv2rgb(yuv):
    """YUV to RGB color space conversion.

    Parameters
    ----------
    yuv : array_like
        The image in YUV format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `yuv` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YUV

    """
    return _convert(rgb_from_yuv, yuv)


def yiq2rgb(yiq):
    """YIQ to RGB color space conversion.

    Parameters
    ----------
    yiq : array_like
        The image in YIQ format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `yiq` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.
    """
    return _convert(rgb_from_yiq, yiq)


def ypbpr2rgb(ypbpr):
    """YPbPr to RGB color space conversion.

    Parameters
    ----------
    ypbpr : array_like
        The image in YPbPr format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `ypbpr` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YPbPr

    """
    return _convert(rgb_from_ypbpr, ypbpr)


def ycbcr2rgb(ycbcr):
    """YCbCr to RGB color space conversion.

    Parameters
    ----------
    ycbcr : array_like
        The image in YCbCr format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `ycbcr` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    Notes
    -----
    Y is between 16 and 235.  This is the color space which is commonly used
    by video codecs, it is sometimes incorrectly called "YUV"

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YCbCr

    """
    arr = ycbcr.copy()
    arr[..., 0] -= 16
    arr[..., 1] -= 128
    arr[..., 2] -= 128
    return _convert(rgb_from_ycbcr, arr)


def ydbdr2rgb(ydbdr):
    """YDbDr to RGB color space conversion.

    Parameters
    ----------
    ydbdr : array_like
        The image in YDbDr format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Returns
    -------
    out : ndarray
        The image in RGB format, in a 3- or 4-D array of shape
        ``(M, N, [P,] 3)``.

    Raises
    ------
    ValueError
        If `ydbdr` is not a 3- or 4-D array of shape ``(M, N, [P,] 3)``.

    Notes
    -----
    This is the color space which is commonly used
    by video codecs, it is also the reversible color transform in JPEG2000.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/YDbDr

    """
    arr = ydbdr.copy()
    return _convert(rgb_from_ydbdr, arr)
