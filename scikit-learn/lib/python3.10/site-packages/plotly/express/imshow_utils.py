"""Vendored code from scikit-image in order to limit the number of dependencies
Extracted from scikit-image/skimage/exposure/exposure.py
"""

import numpy as np

from warnings import warn

_integer_types = (
    np.byte,
    np.ubyte,  # 8 bits
    np.short,
    np.ushort,  # 16 bits
    np.intc,
    np.uintc,  # 16 or 32 or 64 bits
    np.int_,
    np.uint,  # 32 or 64 bits
    np.longlong,
    np.ulonglong,
)  # 64 bits
_integer_ranges = {t: (np.iinfo(t).min, np.iinfo(t).max) for t in _integer_types}
dtype_range = {
    np.bool_: (False, True),
    np.float16: (-1, 1),
    np.float32: (-1, 1),
    np.float64: (-1, 1),
}
dtype_range.update(_integer_ranges)


DTYPE_RANGE = dtype_range.copy()
DTYPE_RANGE.update((d.__name__, limits) for d, limits in dtype_range.items())
DTYPE_RANGE.update(
    {
        "uint10": (0, 2**10 - 1),
        "uint12": (0, 2**12 - 1),
        "uint14": (0, 2**14 - 1),
        "bool": dtype_range[np.bool_],
        "float": dtype_range[np.float64],
    }
)


def intensity_range(image, range_values="image", clip_negative=False):
    """Return image intensity range (min, max) based on desired value type.

    Parameters
    ----------
    image : array
        Input image.
    range_values : str or 2-tuple, optional
        The image intensity range is configured by this parameter.
        The possible values for this parameter are enumerated below.

        'image'
            Return image min/max as the range.
        'dtype'
            Return min/max of the image's dtype as the range.
        dtype-name
            Return intensity range based on desired `dtype`. Must be valid key
            in `DTYPE_RANGE`. Note: `image` is ignored for this range type.
        2-tuple
            Return `range_values` as min/max intensities. Note that there's no
            reason to use this function if you just want to specify the
            intensity range explicitly. This option is included for functions
            that use `intensity_range` to support all desired range types.

    clip_negative : bool, optional
        If True, clip the negative range (i.e. return 0 for min intensity)
        even if the image dtype allows negative values.
    """
    if range_values == "dtype":
        range_values = image.dtype.type

    if range_values == "image":
        i_min = np.min(image)
        i_max = np.max(image)
    elif range_values in DTYPE_RANGE:
        i_min, i_max = DTYPE_RANGE[range_values]
        if clip_negative:
            i_min = 0
    else:
        i_min, i_max = range_values
    return i_min, i_max


def _output_dtype(dtype_or_range):
    """Determine the output dtype for rescale_intensity.

    The dtype is determined according to the following rules:
    - if ``dtype_or_range`` is a dtype, that is the output dtype.
    - if ``dtype_or_range`` is a dtype string, that is the dtype used, unless
      it is not a NumPy data type (e.g. 'uint12' for 12-bit unsigned integers),
      in which case the data type that can contain it will be used
      (e.g. uint16 in this case).
    - if ``dtype_or_range`` is a pair of values, the output data type will be
      float.

    Parameters
    ----------
    dtype_or_range : type, string, or 2-tuple of int/float
        The desired range for the output, expressed as either a NumPy dtype or
        as a (min, max) pair of numbers.

    Returns
    -------
    out_dtype : type
        The data type appropriate for the desired output.
    """
    if type(dtype_or_range) in [list, tuple, np.ndarray]:
        # pair of values: always return float.
        return np.float_
    if isinstance(dtype_or_range, type):
        # already a type: return it
        return dtype_or_range
    if dtype_or_range in DTYPE_RANGE:
        # string key in DTYPE_RANGE dictionary
        try:
            # if it's a canonical numpy dtype, convert
            return np.dtype(dtype_or_range).type
        except TypeError:  # uint10, uint12, uint14
            # otherwise, return uint16
            return np.uint16
    else:
        raise ValueError(
            "Incorrect value for out_range, should be a valid image data "
            "type or a pair of values, got %s." % str(dtype_or_range)
        )


def rescale_intensity(image, in_range="image", out_range="dtype"):
    """Return image after stretching or shrinking its intensity levels.

    The desired intensity range of the input and output, `in_range` and
    `out_range` respectively, are used to stretch or shrink the intensity range
    of the input image. See examples below.

    Parameters
    ----------
    image : array
        Image array.
    in_range, out_range : str or 2-tuple, optional
        Min and max intensity values of input and output image.
        The possible values for this parameter are enumerated below.

        'image'
            Use image min/max as the intensity range.
        'dtype'
            Use min/max of the image's dtype as the intensity range.
        dtype-name
            Use intensity range based on desired `dtype`. Must be valid key
            in `DTYPE_RANGE`.
        2-tuple
            Use `range_values` as explicit min/max intensities.

    Returns
    -------
    out : array
        Image array after rescaling its intensity. This image is the same dtype
        as the input image.

    Notes
    -----
    .. versionchanged:: 0.17
        The dtype of the output array has changed to match the output dtype, or
        float if the output range is specified by a pair of floats.

    See Also
    --------
    equalize_hist

    Examples
    --------
    By default, the min/max intensities of the input image are stretched to
    the limits allowed by the image's dtype, since `in_range` defaults to
    'image' and `out_range` defaults to 'dtype':

    >>> image = np.array([51, 102, 153], dtype=np.uint8)
    >>> rescale_intensity(image)
    array([  0, 127, 255], dtype=uint8)

    It's easy to accidentally convert an image dtype from uint8 to float:

    >>> 1.0 * image
    array([ 51., 102., 153.])

    Use `rescale_intensity` to rescale to the proper range for float dtypes:

    >>> image_float = 1.0 * image
    >>> rescale_intensity(image_float)
    array([0. , 0.5, 1. ])

    To maintain the low contrast of the original, use the `in_range` parameter:

    >>> rescale_intensity(image_float, in_range=(0, 255))
    array([0.2, 0.4, 0.6])

    If the min/max value of `in_range` is more/less than the min/max image
    intensity, then the intensity levels are clipped:

    >>> rescale_intensity(image_float, in_range=(0, 102))
    array([0.5, 1. , 1. ])

    If you have an image with signed integers but want to rescale the image to
    just the positive range, use the `out_range` parameter. In that case, the
    output dtype will be float:

    >>> image = np.array([-10, 0, 10], dtype=np.int8)
    >>> rescale_intensity(image, out_range=(0, 127))
    array([  0. ,  63.5, 127. ])

    To get the desired range with a specific dtype, use ``.astype()``:

    >>> rescale_intensity(image, out_range=(0, 127)).astype(np.int8)
    array([  0,  63, 127], dtype=int8)

    If the input image is constant, the output will be clipped directly to the
    output range:
    >>> image = np.array([130, 130, 130], dtype=np.int32)
    >>> rescale_intensity(image, out_range=(0, 127)).astype(np.int32)
    array([127, 127, 127], dtype=int32)
    """
    if out_range in ["dtype", "image"]:
        out_dtype = _output_dtype(image.dtype.type)
    else:
        out_dtype = _output_dtype(out_range)

    imin, imax = map(float, intensity_range(image, in_range))
    omin, omax = map(
        float, intensity_range(image, out_range, clip_negative=(imin >= 0))
    )

    if np.any(np.isnan([imin, imax, omin, omax])):
        warn(
            "One or more intensity levels are NaN. Rescaling will broadcast "
            "NaN to the full image. Provide intensity levels yourself to "
            "avoid this. E.g. with np.nanmin(image), np.nanmax(image).",
            stacklevel=2,
        )

    image = np.clip(image, imin, imax)

    if imin != imax:
        image = (image - imin) / (imax - imin)
        return np.asarray(image * (omax - omin) + omin, dtype=out_dtype)
    else:
        return np.clip(image, omin, omax).astype(out_dtype)
