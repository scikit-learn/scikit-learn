"""
Utility functions wrapping PIL functions

This is a local version of utility functions from scipy that are wrapping PIL
functionality. These functions are deprecated in scipy 1.0.0 and will be
removed in scipy 1.2.0. Therefore, the functionality used in sklearn is
copied here. Origin is the file scipy/misc/pilutil.py. Parameters and
functions that are not used in sklearn were removed.

Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2017 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Enthought nor the names of the SciPy Developers
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
"""
from __future__ import division, print_function, absolute_import


__all__ = ['_have_image', '_bytescale', '_imread', '_imsave',
           '_fromimage', '_toimage', '_imresize']


import numpy as np

_have_image = True
try:
    try:
        from PIL import Image
    except ImportError:
        import Image
    if not hasattr(Image, 'frombytes'):
        Image.frombytes = Image.fromstring
except ImportError:
    _have_image = False


def _bytescale(data):
    """
    Byte scales an array (image).

    Byte scaling means converting the input image to uint8 dtype and scaling
    the range to ``(0, 255)``.

    If the input image already has dtype uint8, no scaling is done.
    This function is only available if Python Imaging Library (PIL) is
    installed.

    Parameters
    ----------
    data : ndarray
        PIL image data array.

    Returns
    -------
    img_array : uint8 ndarray
        The byte-scaled array.

    Examples
    --------
    >>> from sklearn.externals.pilutil import _bytescale
    >>> img = np.array([[ 91.06794177,   3.39058326,  84.4221549 ],
    ...                 [ 73.88003259,  80.91433048,   4.88878881],
    ...                 [ 51.53875334,  34.45808177,  27.5873488 ]])
    >>> _bytescale(img)
    array([[255,   0, 236],
           [205, 225,   4],
           [140,  90,  70]], dtype=uint8)
    """
    if data.dtype == np.uint8:
        return data

    cmin = data.min()
    cmax = data.max()

    cscale = cmax - cmin
    if cscale == 0:
        cscale = 1

    scale = 255. / cscale
    bytedata = (data - cmin) * scale
    return (bytedata.clip(0, 255) + 0.5).astype(np.uint8)


def _imread(name):
    """
    Read an image from a file as an array.

    This function is only available if Python Imaging Library (PIL) is
    installed.

    Parameters
    ----------
    name : str or file object
        The file name or file object to be read.

    Returns
    -------
    imread : ndarray
        The array obtained by reading the image.

    Notes
    -----
    This is a simplified combination of scipy's scipy.misc.pilutil.imread and
    scipy.misc.pilutil.fromimage, which are deprecated in scipy 1.0.0 and will
    be removed from scipy in version 1.2.0.
    """
    if not _have_image:
        raise ImportError("The Python Imaging Library (PIL) "
                          "is required to load data from jpeg files")
    pil_image = Image.open(name)

    return _fromimage(pil_image)


def _imsave(name, arr):
    """
    Save an array as an image.
    This function is only available if Python Imaging Library (PIL) is installed.
    .. warning::
        This function uses `bytescale` under the hood to rescale images to use
        the full (0, 255) range if ``mode`` is one of ``None, 'L', 'P', 'l'``.
        It will also cast data for 2-D images to ``uint32`` for ``mode=None``
        (which is the default).

    Parameters
    ----------
    name : str or file object
        Output file name or file object.
    arr : ndarray, MxN or MxNx3 or MxNx4
        Array containing image values.  If the shape is ``MxN``, the array
        represents a grey-level image.  Shape ``MxNx3`` stores the red, green
        and blue bands along the last dimension.  An alpha layer may be
        included, specified as the last colour band of an ``MxNx4`` array.
    """
    pil_image = _toimage(arr, channel_axis=2)
    pil_image.save(name)
    return


def _fromimage(pil_image):
    """
    Return a copy of a PIL image as a numpy array.

    This function is only available if Python Imaging Library (PIL) is
    installed.

    Parameters
    ----------
    im : PIL image
        Input image.

    Returns
    -------
    fromimage : ndarray
        The different colour bands/channels are stored in the
        third dimension, such that a grey-image is MxN, an
        RGB-image MxNx3 and an RGBA-image MxNx4.
    """
    if not _have_image:
        raise ImportError("The Python Imaging Library (PIL) "
                          "is required to load data from jpeg files")
    if not Image.isImageType(pil_image):
        raise TypeError("Input is not a PIL image.")

    if pil_image.mode == 'P':
        # Mode 'P' means there is an indexed "palette".  If we leave the mode
        # as 'P', then when we do `a = array(pil_image)` below, `a` will be a
        # 2-D containing the indices into the palette, and not a 3-D array
        # containing the RGB or RGBA values.
        if 'transparency' in pil_image.info:
            pil_image = pil_image.convert('RGBA')
        else:
            pil_image = pil_image.convert('RGB')

    if pil_image.mode == '1':
        # Workaround for crash in PIL. When pil_image is 1-bit, the cal
        # array(pil_image) can cause a seg. fault, or generate garbage. See
        # https://github.com/scipy/scipy/issues/2138 and
        # https://github.com/python-pillow/Pillow/issues/350.
        #
        # This converts im from a 1-bit image to an 8-bit image.
        pil_image = pil_image.convert('L')

    return np.array(pil_image)


def _toimage(arr, channel_axis=None):
    """
    Takes a numpy array and returns a PIL image.

    This function is only available if Python Imaging Library (PIL) is
    installed.
    .. warning::
        This function uses `_bytescale` under the hood to rescale images to
        use the full (0, 255) range. It will also cast data for 2-D images to
        ``uint32``.

    Notes
    -----
    For 3-D arrays if one of the dimensions is 3, the mode is 'RGB'
    by default or 'YCbCr' if selected.
    The numpy array must be either 2 dimensional or 3 dimensional.
    """
    if not _have_image:
        raise ImportError("The Python Imaging Library (PIL) "
                          "is required to load data from jpeg files")
    data = np.asarray(arr)
    if np.iscomplexobj(data):
        raise ValueError("Cannot convert a complex-valued array.")
    shape = list(data.shape)
    valid = len(shape) == 2 or ((len(shape) == 3) and
                                ((3 in shape) or (4 in shape)))
    if not valid:
        raise ValueError("'arr' does not have a suitable array shape for "
                         "any mode.")
    if len(shape) == 2:
        shape = (shape[1], shape[0])
        bytedata = _bytescale(data)
        image = Image.frombytes('L', shape, bytedata.tostring())
        return image

    # if here then 3-d array with a 3 or a 4 in the shape length.
    # Check for 3 in datacube shape --- 'RGB' or 'YCbCr'
    if channel_axis is None:
        if 3 in shape:
            ca = np.flatnonzero(np.asarray(shape) == 3)[0]
        else:
            ca = np.flatnonzero(np.asarray(shape) == 4)
            if not ca:
                ca = ca[0]
            else:
                raise ValueError("Could not find channel dimension.")
    else:
        ca = channel_axis

    numch = shape[ca]
    if numch not in [3, 4]:
        raise ValueError("Channel axis dimension is not valid.")

    bytedata = _bytescale(data)
    if ca == 2:
        strdata = bytedata.tostring()
        shape = (shape[1], shape[0])
    elif ca == 1:
        strdata = np.transpose(bytedata, (0, 2, 1)).tostring()
        shape = (shape[2], shape[0])
    elif ca == 0:
        strdata = np.transpose(bytedata, (1, 2, 0)).tostring()
        shape = (shape[2], shape[1])
    else:
        raise ValueError("Invalid channel dimension.")

    if numch == 3:
        mode = 'RGB'
    else:
        mode = 'RGBA'

    # Here we know data and mode is correct
    return Image.frombytes(mode, shape, strdata)


def _imresize(arr, size):
    """
    Resize an image.

    This function is only available if Python Imaging Library (PIL) is
    installed.
    .. warning::
        This function uses `_bytescale` under the hood to rescale images to
        use the full (0, 255) range.
        It will also cast data for 2-D images to ``uint32``.

    Parameters
    ----------
    arr : ndarray
        The array of image to be resized.
    size : int, float or tuple
        * int   - Percentage of current size.
        * float - Fraction of current size.
        * tuple - Size of the output image (height, width).

    Returns
    -------
    imresize : ndarray
        The resized array of image.
    """
    im = _toimage(arr)
    ts = type(size)
    if np.issubdtype(ts, np.signedinteger):
        percent = size / 100.0
        size = tuple((np.array(im.size) * percent).astype(int))
    elif np.issubdtype(type(size), np.floating):
        size = tuple((np.array(im.size) * size).astype(int))
    else:
        size = (size[1], size[0])
    imnew = im.resize(size, resample=2)
    return _fromimage(imnew)
