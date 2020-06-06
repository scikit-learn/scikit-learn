"""Docstring components common to several ndimage functions."""
from __future__ import division, print_function, absolute_import

from scipy._lib import doccer

__all__ = ['docfiller']


_input_doc = (
"""input : array_like
    The input array.""")
_axis_doc = (
"""axis : int, optional
    The axis of `input` along which to calculate. Default is -1.""")
_output_doc = (
"""output : array or dtype, optional
    The array in which to place the output, or the dtype of the
    returned array. By default an array of the same dtype as input
    will be created.""")
_size_foot_doc = (
"""size : scalar or tuple, optional
    See footprint, below. Ignored if footprint is given.
footprint : array, optional
    Either `size` or `footprint` must be defined.  `size` gives
    the shape that is taken from the input array, at every element
    position, to define the input to the filter function.
    `footprint` is a boolean array that specifies (implicitly) a
    shape, but also which of the elements within this shape will get
    passed to the filter function.  Thus ``size=(n,m)`` is equivalent
    to ``footprint=np.ones((n,m))``.  We adjust `size` to the number
    of dimensions of the input array, so that, if the input array is
    shape (10,10,10), and `size` is 2, then the actual size used is
    (2,2,2). When `footprint` is given, `size` is ignored.""")
_mode_doc = (
"""mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
    The `mode` parameter determines how the input array is extended
    beyond its boundaries. Default is 'reflect'. Behavior for each valid
    value is as follows:

    'reflect' (`d c b a | a b c d | d c b a`)
        The input is extended by reflecting about the edge of the last
        pixel.

    'constant' (`k k k k | a b c d | k k k k`)
        The input is extended by filling all values beyond the edge with
        the same constant value, defined by the `cval` parameter.

    'nearest' (`a a a a | a b c d | d d d d`)
        The input is extended by replicating the last pixel.

    'mirror' (`d c b | a b c d | c b a`)
        The input is extended by reflecting about the center of the last
        pixel.

    'wrap' (`a b c d | a b c d | a b c d`)
        The input is extended by wrapping around to the opposite edge.""")
_mode_multiple_doc = (
"""mode : str or sequence, optional
    The `mode` parameter determines how the input array is extended
    when the filter overlaps a border. By passing a sequence of modes
    with length equal to the number of dimensions of the input array,
    different modes can be specified along each axis. Default value is
    'reflect'. The valid values and their behavior is as follows:

    'reflect' (`d c b a | a b c d | d c b a`)
        The input is extended by reflecting about the edge of the last
        pixel.

    'constant' (`k k k k | a b c d | k k k k`)
        The input is extended by filling all values beyond the edge with
        the same constant value, defined by the `cval` parameter.

    'nearest' (`a a a a | a b c d | d d d d`)
        The input is extended by replicating the last pixel.

    'mirror' (`d c b | a b c d | c b a`)
        The input is extended by reflecting about the center of the last
        pixel.

    'wrap' (`a b c d | a b c d | a b c d`)
        The input is extended by wrapping around to the opposite edge.""")
_cval_doc = (
"""cval : scalar, optional
    Value to fill past edges of input if `mode` is 'constant'. Default
    is 0.0.""")
_origin_doc = (
"""origin : int, optional
    Controls the placement of the filter on the input array's pixels.
    A value of 0 (the default) centers the filter over the pixel, with
    positive values shifting the filter to the left, and negative ones
    to the right.""")
_origin_multiple_doc = (
"""origin : int or sequence, optional
    Controls the placement of the filter on the input array's pixels.
    A value of 0 (the default) centers the filter over the pixel, with
    positive values shifting the filter to the left, and negative ones
    to the right. By passing a sequence of origins with length equal to
    the number of dimensions of the input array, different shifts can
    be specified along each axis.""")
_extra_arguments_doc = (
"""extra_arguments : sequence, optional
    Sequence of extra positional arguments to pass to passed function.""")
_extra_keywords_doc = (
"""extra_keywords : dict, optional
    dict of extra keyword arguments to pass to passed function.""")
_prefilter_doc = (
"""prefilter : bool, optional
    Determines if the input array is prefiltered with `spline_filter`
    before interpolation. The default is True, which will create a
    temporary `float64` array of filtered values if `order > 1`. If
    setting this to False, the output will be slightly blurred if
    `order > 1`, unless the input is prefiltered, i.e. it is the result
    of calling `spline_filter` on the original input.""")

docdict = {
    'input': _input_doc,
    'axis': _axis_doc,
    'output': _output_doc,
    'size_foot': _size_foot_doc,
    'mode': _mode_doc,
    'mode_multiple': _mode_multiple_doc,
    'cval': _cval_doc,
    'origin': _origin_doc,
    'origin_multiple': _origin_multiple_doc,
    'extra_arguments': _extra_arguments_doc,
    'extra_keywords': _extra_keywords_doc,
    'prefilter': _prefilter_doc
    }

docfiller = doccer.filldoc(docdict)
