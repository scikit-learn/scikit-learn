#
# The Python Imaging Library.
# $Id$
#
# standard filters
#
# History:
# 1995-11-27 fl   Created
# 2002-06-08 fl   Added rank and mode filters
# 2003-09-15 fl   Fixed rank calculation in rank filter; added expand call
#
# Copyright (c) 1997-2003 by Secret Labs AB.
# Copyright (c) 1995-2002 by Fredrik Lundh.
#
# See the README file for information on usage and redistribution.
#

import functools


class Filter(object):
    pass


class MultibandFilter(Filter):
    pass


class Kernel(MultibandFilter):
    """
    Create a convolution kernel.  The current version only
    supports 3x3 and 5x5 integer and floating point kernels.

    In the current version, kernels can only be applied to
    "L" and "RGB" images.

    :param size: Kernel size, given as (width, height). In the current
                    version, this must be (3,3) or (5,5).
    :param kernel: A sequence containing kernel weights.
    :param scale: Scale factor. If given, the result for each pixel is
                    divided by this value.  the default is the sum of the
                    kernel weights.
    :param offset: Offset. If given, this value is added to the result,
                    after it has been divided by the scale factor.
    """

    def __init__(self, size, kernel, scale=None, offset=0):
        if scale is None:
            # default scale is sum of kernel
            scale = functools.reduce(lambda a, b: a+b, kernel)
        if size[0] * size[1] != len(kernel):
            raise ValueError("not enough coefficients in kernel")
        self.filterargs = size, scale, offset, kernel

    def filter(self, image):
        if image.mode == "P":
            raise ValueError("cannot filter palette images")
        return image.filter(*self.filterargs)


class BuiltinFilter(Kernel):
    def __init__(self):
        pass


class RankFilter(Filter):
    """
    Create a rank filter.  The rank filter sorts all pixels in
    a window of the given size, and returns the **rank**'th value.

    :param size: The kernel size, in pixels.
    :param rank: What pixel value to pick.  Use 0 for a min filter,
                 ``size * size / 2`` for a median filter, ``size * size - 1``
                 for a max filter, etc.
    """
    name = "Rank"

    def __init__(self, size, rank):
        self.size = size
        self.rank = rank

    def filter(self, image):
        if image.mode == "P":
            raise ValueError("cannot filter palette images")
        image = image.expand(self.size//2, self.size//2)
        return image.rankfilter(self.size, self.rank)


class MedianFilter(RankFilter):
    """
    Create a median filter. Picks the median pixel value in a window with the
    given size.

    :param size: The kernel size, in pixels.
    """
    name = "Median"

    def __init__(self, size=3):
        self.size = size
        self.rank = size*size//2


class MinFilter(RankFilter):
    """
    Create a min filter.  Picks the lowest pixel value in a window with the
    given size.

    :param size: The kernel size, in pixels.
    """
    name = "Min"

    def __init__(self, size=3):
        self.size = size
        self.rank = 0


class MaxFilter(RankFilter):
    """
    Create a max filter.  Picks the largest pixel value in a window with the
    given size.

    :param size: The kernel size, in pixels.
    """
    name = "Max"

    def __init__(self, size=3):
        self.size = size
        self.rank = size*size-1


class ModeFilter(Filter):
    """

    Create a mode filter. Picks the most frequent pixel value in a box with the
    given size.  Pixel values that occur only once or twice are ignored; if no
    pixel value occurs more than twice, the original pixel value is preserved.

    :param size: The kernel size, in pixels.
    """
    name = "Mode"

    def __init__(self, size=3):
        self.size = size

    def filter(self, image):
        return image.modefilter(self.size)


class GaussianBlur(MultibandFilter):
    """Gaussian blur filter.

    :param radius: Blur radius.
    """
    name = "GaussianBlur"

    def __init__(self, radius=2):
        self.radius = radius

    def filter(self, image):
        return image.gaussian_blur(self.radius)


class BoxBlur(MultibandFilter):
    """Blurs the image by setting each pixel to the average value of the pixels
    in a square box extending radius pixels in each direction.
    Supports float radius of arbitrary size. Uses an optimized implementation
    which runs in linear time relative to the size of the image
    for any radius value.

    :param radius: Size of the box in one direction. Radius 0 does not blur,
                   returns an identical image. Radius 1 takes 1 pixel
                   in each direction, i.e. 9 pixels in total.
    """
    name = "BoxBlur"

    def __init__(self, radius):
        self.radius = radius

    def filter(self, image):
        return image.box_blur(self.radius)


class UnsharpMask(MultibandFilter):
    """Unsharp mask filter.

    See Wikipedia's entry on `digital unsharp masking`_ for an explanation of
    the parameters.

    :param radius: Blur Radius
    :param percent: Unsharp strength, in percent
    :param threshold: Threshold controls the minimum brightness change that
      will be sharpened

    .. _digital unsharp masking: https://en.wikipedia.org/wiki/Unsharp_masking#Digital_unsharp_masking

    """
    name = "UnsharpMask"

    def __init__(self, radius=2, percent=150, threshold=3):
        self.radius = radius
        self.percent = percent
        self.threshold = threshold

    def filter(self, image):
        return image.unsharp_mask(self.radius, self.percent, self.threshold)


class BLUR(BuiltinFilter):
    name = "Blur"
    filterargs = (5, 5), 16, 0, (
        1,  1,  1,  1,  1,
        1,  0,  0,  0,  1,
        1,  0,  0,  0,  1,
        1,  0,  0,  0,  1,
        1,  1,  1,  1,  1
        )


class CONTOUR(BuiltinFilter):
    name = "Contour"
    filterargs = (3, 3), 1, 255, (
        -1, -1, -1,
        -1,  8, -1,
        -1, -1, -1
        )


class DETAIL(BuiltinFilter):
    name = "Detail"
    filterargs = (3, 3), 6, 0, (
        0, -1,  0,
        -1, 10, -1,
        0, -1,  0
        )


class EDGE_ENHANCE(BuiltinFilter):
    name = "Edge-enhance"
    filterargs = (3, 3), 2, 0, (
        -1, -1, -1,
        -1, 10, -1,
        -1, -1, -1
        )


class EDGE_ENHANCE_MORE(BuiltinFilter):
    name = "Edge-enhance More"
    filterargs = (3, 3), 1, 0, (
        -1, -1, -1,
        -1,  9, -1,
        -1, -1, -1
        )


class EMBOSS(BuiltinFilter):
    name = "Emboss"
    filterargs = (3, 3), 1, 128, (
        -1,  0,  0,
        0,  1,  0,
        0,  0,  0
        )


class FIND_EDGES(BuiltinFilter):
    name = "Find Edges"
    filterargs = (3, 3), 1, 0, (
        -1, -1, -1,
        -1,  8, -1,
        -1, -1, -1
        )


class SHARPEN(BuiltinFilter):
    name = "Sharpen"
    filterargs = (3, 3), 16, 0, (
        -2, -2, -2,
        -2, 32, -2,
        -2, -2, -2
        )


class SMOOTH(BuiltinFilter):
    name = "Smooth"
    filterargs = (3, 3), 13, 0, (
        1,  1,  1,
        1,  5,  1,
        1,  1,  1
        )


class SMOOTH_MORE(BuiltinFilter):
    name = "Smooth More"
    filterargs = (5, 5), 100, 0, (
        1,  1,  1,  1,  1,
        1,  5,  5,  5,  1,
        1,  5, 44,  5,  1,
        1,  5,  5,  5,  1,
        1,  1,  1,  1,  1
        )
