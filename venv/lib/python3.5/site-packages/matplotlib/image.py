"""
The image module supports basic image loading, rescaling and display
operations.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves.urllib.parse import urlparse
from six.moves.urllib.request import urlopen
from io import BytesIO

from math import ceil
import os
import logging

import numpy as np

from matplotlib import rcParams
import matplotlib.artist as martist
from matplotlib.artist import allow_rasterization
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.cbook as cbook

# For clarity, names from _image are given explicitly in this module:
import matplotlib._image as _image
import matplotlib._png as _png

# For user convenience, the names from _image are also imported into
# the image namespace:
from matplotlib._image import *

from matplotlib.transforms import (Affine2D, BboxBase, Bbox, BboxTransform,
                                   IdentityTransform, TransformedBbox)

_log = logging.getLogger(__name__)

# map interpolation strings to module constants
_interpd_ = {
    'none': _image.NEAREST,  # fall back to nearest when not supported
    'nearest': _image.NEAREST,
    'bilinear': _image.BILINEAR,
    'bicubic': _image.BICUBIC,
    'spline16': _image.SPLINE16,
    'spline36': _image.SPLINE36,
    'hanning': _image.HANNING,
    'hamming': _image.HAMMING,
    'hermite': _image.HERMITE,
    'kaiser': _image.KAISER,
    'quadric': _image.QUADRIC,
    'catrom': _image.CATROM,
    'gaussian': _image.GAUSSIAN,
    'bessel': _image.BESSEL,
    'mitchell': _image.MITCHELL,
    'sinc': _image.SINC,
    'lanczos': _image.LANCZOS,
    'blackman': _image.BLACKMAN,
}

interpolations_names = set(_interpd_)


def composite_images(images, renderer, magnification=1.0):
    """
    Composite a number of RGBA images into one.  The images are
    composited in the order in which they appear in the `images` list.

    Parameters
    ----------
    images : list of Images
        Each must have a `make_image` method.  For each image,
        `can_composite` should return `True`, though this is not
        enforced by this function.  Each image must have a purely
        affine transformation with no shear.

    renderer : RendererBase instance

    magnification : float
        The additional magnification to apply for the renderer in use.

    Returns
    -------
    tuple : image, offset_x, offset_y
        Returns the tuple:

        - image: A numpy array of the same type as the input images.

        - offset_x, offset_y: The offset of the image (left, bottom)
          in the output figure.
    """
    if len(images) == 0:
        return np.empty((0, 0, 4), dtype=np.uint8), 0, 0

    parts = []
    bboxes = []
    for image in images:
        data, x, y, trans = image.make_image(renderer, magnification)
        if data is not None:
            x *= magnification
            y *= magnification
            parts.append((data, x, y, image.get_alpha() or 1.0))
            bboxes.append(
                Bbox([[x, y], [x + data.shape[1], y + data.shape[0]]]))

    if len(parts) == 0:
        return np.empty((0, 0, 4), dtype=np.uint8), 0, 0

    bbox = Bbox.union(bboxes)

    output = np.zeros(
        (int(bbox.height), int(bbox.width), 4), dtype=np.uint8)

    for data, x, y, alpha in parts:
        trans = Affine2D().translate(x - bbox.x0, y - bbox.y0)
        _image.resample(data, output, trans, _image.NEAREST,
                        resample=False, alpha=alpha)

    return output, bbox.x0 / magnification, bbox.y0 / magnification


def _draw_list_compositing_images(
        renderer, parent, artists, suppress_composite=None):
    """
    Draw a sorted list of artists, compositing images into a single
    image where possible.

    For internal matplotlib use only: It is here to reduce duplication
    between `Figure.draw` and `Axes.draw`, but otherwise should not be
    generally useful.
    """
    has_images = any(isinstance(x, _ImageBase) for x in artists)

    # override the renderer default if suppressComposite is not None
    not_composite = (suppress_composite if suppress_composite is not None
                     else renderer.option_image_nocomposite())

    if not_composite or not has_images:
        for a in artists:
            a.draw(renderer)
    else:
        # Composite any adjacent images together
        image_group = []
        mag = renderer.get_image_magnification()

        def flush_images():
            if len(image_group) == 1:
                image_group[0].draw(renderer)
            elif len(image_group) > 1:
                data, l, b = composite_images(image_group, renderer, mag)
                if data.size != 0:
                    gc = renderer.new_gc()
                    gc.set_clip_rectangle(parent.bbox)
                    gc.set_clip_path(parent.get_clip_path())
                    renderer.draw_image(gc, np.round(l), np.round(b), data)
                    gc.restore()
            del image_group[:]

        for a in artists:
            if isinstance(a, _ImageBase) and a.can_composite():
                image_group.append(a)
            else:
                flush_images()
                a.draw(renderer)
        flush_images()


def _rgb_to_rgba(A):
    """
    Convert an RGB image to RGBA, as required by the image resample C++
    extension.
    """
    rgba = np.zeros((A.shape[0], A.shape[1], 4), dtype=A.dtype)
    rgba[:, :, :3] = A
    if rgba.dtype == np.uint8:
        rgba[:, :, 3] = 255
    else:
        rgba[:, :, 3] = 1.0
    return rgba


class _ImageBase(martist.Artist, cm.ScalarMappable):
    zorder = 0

    @property
    @cbook.deprecated("2.1")
    def _interpd(self):
        return _interpd_

    @property
    @cbook.deprecated("2.1")
    def _interpdr(self):
        return {v: k for k, v in six.iteritems(_interpd_)}

    @property
    @cbook.deprecated("2.1", alternative="mpl.image.interpolation_names")
    def iterpnames(self):
        return interpolations_names

    def __str__(self):
        return "AxesImage(%g,%g;%gx%g)" % tuple(self.axes.bbox.bounds)

    def __init__(self, ax,
                 cmap=None,
                 norm=None,
                 interpolation=None,
                 origin=None,
                 filternorm=1,
                 filterrad=4.0,
                 resample=False,
                 **kwargs
                 ):
        """
        interpolation and cmap default to their rc settings

        cmap is a colors.Colormap instance
        norm is a colors.Normalize instance to map luminance to 0-1

        extent is data axes (left, right, bottom, top) for making image plots
        registered with data plots.  Default is to label the pixel
        centers with the zero-based row and column indices.

        Additional kwargs are matplotlib.artist properties

        """
        martist.Artist.__init__(self)
        cm.ScalarMappable.__init__(self, norm, cmap)
        self._mouseover = True
        if origin is None:
            origin = rcParams['image.origin']
        self.origin = origin
        self.set_filternorm(filternorm)
        self.set_filterrad(filterrad)
        self.set_interpolation(interpolation)
        self.set_resample(resample)
        self.axes = ax

        self._imcache = None

        self.update(kwargs)

    def __getstate__(self):
        state = super(_ImageBase, self).__getstate__()
        # We can't pickle the C Image cached object.
        state['_imcache'] = None
        return state

    def get_size(self):
        """Get the numrows, numcols of the input image"""
        if self._A is None:
            raise RuntimeError('You must first set the image array')

        return self._A.shape[:2]

    def set_alpha(self, alpha):
        """
        Set the alpha value used for blending - not supported on
        all backends

        ACCEPTS: float
        """
        martist.Artist.set_alpha(self, alpha)
        self._imcache = None

    def changed(self):
        """
        Call this whenever the mappable is changed so observers can
        update state
        """
        self._imcache = None
        self._rgbacache = None
        cm.ScalarMappable.changed(self)

    def _make_image(self, A, in_bbox, out_bbox, clip_bbox, magnification=1.0,
                    unsampled=False, round_to_pixel_border=True):
        """
        Normalize, rescale and color the image `A` from the given
        in_bbox (in data space), to the given out_bbox (in pixel
        space) clipped to the given clip_bbox (also in pixel space),
        and magnified by the magnification factor.

        `A` may be a greyscale image (MxN) with a dtype of `float32`,
        `float64`, `uint16` or `uint8`, or an RGBA image (MxNx4) with
        a dtype of `float32`, `float64`, or `uint8`.

        If `unsampled` is True, the image will not be scaled, but an
        appropriate affine transformation will be returned instead.

        If `round_to_pixel_border` is True, the output image size will
        be rounded to the nearest pixel boundary.  This makes the
        images align correctly with the axes.  It should not be used
        in cases where you want exact scaling, however, such as
        FigureImage.

        Returns the resulting (image, x, y, trans), where (x, y) is
        the upper left corner of the result in pixel space, and
        `trans` is the affine transformation from the image to pixel
        space.
        """
        if A is None:
            raise RuntimeError('You must first set the image '
                               'array or the image attribute')
        if A.size == 0:
            raise RuntimeError("_make_image must get a non-empty image. "
                               "Your Artist's draw method must filter before "
                               "this method is called.")

        clipped_bbox = Bbox.intersection(out_bbox, clip_bbox)

        if clipped_bbox is None:
            return None, 0, 0, None

        out_width_base = clipped_bbox.width * magnification
        out_height_base = clipped_bbox.height * magnification

        if out_width_base == 0 or out_height_base == 0:
            return None, 0, 0, None

        if self.origin == 'upper':
            # Flip the input image using a transform.  This avoids the
            # problem with flipping the array, which results in a copy
            # when it is converted to contiguous in the C wrapper
            t0 = Affine2D().translate(0, -A.shape[0]).scale(1, -1)
        else:
            t0 = IdentityTransform()

        t0 += (
            Affine2D()
            .scale(
                in_bbox.width / A.shape[1],
                in_bbox.height / A.shape[0])
            .translate(in_bbox.x0, in_bbox.y0)
            + self.get_transform())

        t = (t0
             + Affine2D().translate(
                 -clipped_bbox.x0,
                 -clipped_bbox.y0)
             .scale(magnification, magnification))

        # So that the image is aligned with the edge of the axes, we want
        # to round up the output width to the next integer.  This also
        # means scaling the transform just slightly to account for the
        # extra subpixel.
        if (t.is_affine and round_to_pixel_border and
                (out_width_base % 1.0 != 0.0 or out_height_base % 1.0 != 0.0)):
            out_width = int(ceil(out_width_base))
            out_height = int(ceil(out_height_base))
            extra_width = (out_width - out_width_base) / out_width_base
            extra_height = (out_height - out_height_base) / out_height_base
            t += Affine2D().scale(1.0 + extra_width, 1.0 + extra_height)
        else:
            out_width = int(out_width_base)
            out_height = int(out_height_base)

        if not unsampled:
            if A.ndim not in (2, 3):
                raise ValueError("Invalid dimensions, got {}".format(A.shape))

            if A.ndim == 2:
                # if we are a 2D array, then we are running through the
                # norm + colormap transformation.  However, in general the
                # input data is not going to match the size on the screen so we
                # have to resample to the correct number of pixels
                # need to

                # TODO slice input array first
                inp_dtype = A.dtype
                a_min = A.min()
                a_max = A.max()
                # figure out the type we should scale to.  For floats,
                # leave as is.  For integers cast to an appropriate-sized
                # float.  Small integers get smaller floats in an attempt
                # to keep the memory footprint reasonable.
                if a_min is np.ma.masked:
                    # all masked, so values don't matter
                    a_min, a_max = np.int32(0), np.int32(1)
                if inp_dtype.kind == 'f':
                    scaled_dtype = A.dtype
                else:
                    # probably an integer of some type.
                    da = a_max.astype(np.float64) - a_min.astype(np.float64)
                    if da > 1e8:
                        # give more breathing room if a big dynamic range
                        scaled_dtype = np.float64
                    else:
                        scaled_dtype = np.float32

                # scale the input data to [.1, .9].  The Agg
                # interpolators clip to [0, 1] internally, use a
                # smaller input scale to identify which of the
                # interpolated points need to be should be flagged as
                # over / under.
                # This may introduce numeric instabilities in very broadly
                # scaled data
                A_scaled = np.empty(A.shape, dtype=scaled_dtype)
                A_scaled[:] = A
                # clip scaled data around norm if necessary.
                # This is necessary for big numbers at the edge of
                # float64's ability to represent changes.  Applying
                # a norm first would be good, but ruins the interpolation
                # of over numbers.
                if self.norm.vmin is not None and self.norm.vmax is not None:
                    dv = (np.float64(self.norm.vmax) -
                            np.float64(self.norm.vmin))
                    vmid = self.norm.vmin + dv / 2
                    newmin = vmid - dv * 1.e7
                    if newmin < a_min:
                        newmin = None
                    else:
                        a_min = np.float64(newmin)
                    newmax = vmid + dv * 1.e7
                    if newmax > a_max:
                        newmax = None
                    else:
                        a_max = np.float64(newmax)
                    if newmax is not None or newmin is not None:
                        A_scaled = np.clip(A_scaled, newmin, newmax)

                A_scaled -= a_min
                # a_min and a_max might be ndarray subclasses so use
                # asscalar to ensure they are scalars to avoid errors
                a_min = np.asscalar(a_min.astype(scaled_dtype))
                a_max = np.asscalar(a_max.astype(scaled_dtype))

                if a_min != a_max:
                    A_scaled /= ((a_max - a_min) / 0.8)
                A_scaled += 0.1
                A_resampled = np.zeros((out_height, out_width),
                                       dtype=A_scaled.dtype)
                # resample the input data to the correct resolution and shape
                _image.resample(A_scaled, A_resampled,
                                t,
                                _interpd_[self.get_interpolation()],
                                self.get_resample(), 1.0,
                                self.get_filternorm() or 0.0,
                                self.get_filterrad() or 0.0)

                # we are done with A_scaled now, remove from namespace
                # to be sure!
                del A_scaled
                # un-scale the resampled data to approximately the
                # original range things that interpolated to above /
                # below the original min/max will still be above /
                # below, but possibly clipped in the case of higher order
                # interpolation + drastically changing data.
                A_resampled -= 0.1
                if a_min != a_max:
                    A_resampled *= ((a_max - a_min) / 0.8)
                A_resampled += a_min
                # if using NoNorm, cast back to the original datatype
                if isinstance(self.norm, mcolors.NoNorm):
                    A_resampled = A_resampled.astype(A.dtype)

                mask = np.empty(A.shape, dtype=np.float32)
                if A.mask.shape == A.shape:
                    # this is the case of a nontrivial mask
                    mask[:] = np.where(A.mask, np.float32(np.nan),
                                       np.float32(1))
                else:
                    mask[:] = 1

                # we always have to interpolate the mask to account for
                # non-affine transformations
                out_mask = np.zeros((out_height, out_width),
                                    dtype=mask.dtype)
                _image.resample(mask, out_mask,
                                t,
                                _interpd_[self.get_interpolation()],
                                True, 1,
                                self.get_filternorm() or 0.0,
                                self.get_filterrad() or 0.0)
                # we are done with the mask, delete from namespace to be sure!
                del mask
                # Agg updates the out_mask in place.  If the pixel has
                # no image data it will not be updated (and still be 0
                # as we initialized it), if input data that would go
                # into that output pixel than it will be `nan`, if all
                # the input data for a pixel is good it will be 1, and
                # if there is _some_ good data in that output pixel it
                # will be between [0, 1] (such as a rotated image).

                out_alpha = np.array(out_mask)
                out_mask = np.isnan(out_mask)
                out_alpha[out_mask] = 1

                # mask and run through the norm
                output = self.norm(np.ma.masked_array(A_resampled, out_mask))
            else:
                # Always convert to RGBA, even if only RGB input
                if A.shape[2] == 3:
                    A = _rgb_to_rgba(A)
                elif A.shape[2] != 4:
                    raise ValueError("Invalid dimensions, got %s" % (A.shape,))

                output = np.zeros((out_height, out_width, 4), dtype=A.dtype)

                alpha = self.get_alpha()
                if alpha is None:
                    alpha = 1.0

                _image.resample(
                    A, output, t, _interpd_[self.get_interpolation()],
                    self.get_resample(), alpha,
                    self.get_filternorm() or 0.0, self.get_filterrad() or 0.0)

            # at this point output is either a 2D array of normed data
            # (of int or float)
            # or an RGBA array of re-sampled input
            output = self.to_rgba(output, bytes=True, norm=False)
            # output is now a correctly sized RGBA array of uint8

            # Apply alpha *after* if the input was greyscale without a mask
            if A.ndim == 2:
                alpha = self.get_alpha()
                if alpha is None:
                    alpha = 1
                alpha_channel = output[:, :, 3]
                alpha_channel[:] = np.asarray(
                    np.asarray(alpha_channel, np.float32) * out_alpha * alpha,
                    np.uint8)

        else:
            if self._imcache is None:
                self._imcache = self.to_rgba(A, bytes=True, norm=(A.ndim == 2))
            output = self._imcache

            # Subset the input image to only the part that will be
            # displayed
            subset = TransformedBbox(
                clip_bbox, t0.frozen().inverted()).frozen()
            output = output[
                int(max(subset.ymin, 0)):
                int(min(subset.ymax + 1, output.shape[0])),
                int(max(subset.xmin, 0)):
                int(min(subset.xmax + 1, output.shape[1]))]

            t = Affine2D().translate(
                int(max(subset.xmin, 0)), int(max(subset.ymin, 0))) + t

        return output, clipped_bbox.x0, clipped_bbox.y0, t

    def make_image(self, renderer, magnification=1.0, unsampled=False):
        raise RuntimeError('The make_image method must be overridden.')

    def _draw_unsampled_image(self, renderer, gc):
        """
        draw unsampled image. The renderer should support a draw_image method
        with scale parameter.
        """

        im, l, b, trans = self.make_image(renderer, unsampled=True)

        if im is None:
            return

        trans = Affine2D().scale(im.shape[1], im.shape[0]) + trans

        renderer.draw_image(gc, l, b, im, trans)

    def _check_unsampled_image(self, renderer):
        """
        return True if the image is better to be drawn unsampled.
        The derived class needs to override it.
        """
        return False

    @allow_rasterization
    def draw(self, renderer, *args, **kwargs):
        # if not visible, declare victory and return
        if not self.get_visible():
            self.stale = False
            return

        # for empty images, there is nothing to draw!
        if self.get_array().size == 0:
            self.stale = False
            return

        # actually render the image.
        gc = renderer.new_gc()
        self._set_gc_clip(gc)
        gc.set_alpha(self.get_alpha())
        gc.set_url(self.get_url())
        gc.set_gid(self.get_gid())

        if (self._check_unsampled_image(renderer) and
                self.get_transform().is_affine):
            self._draw_unsampled_image(renderer, gc)
        else:
            im, l, b, trans = self.make_image(
                renderer, renderer.get_image_magnification())
            if im is not None:
                renderer.draw_image(gc, l, b, im)
        gc.restore()
        self.stale = False

    def contains(self, mouseevent):
        """
        Test whether the mouse event occurred within the image.
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)
        # TODO: make sure this is consistent with patch and patch
        # collection on nonlinear transformed coordinates.
        # TODO: consider returning image coordinates (shouldn't
        # be too difficult given that the image is rectilinear
        x, y = mouseevent.xdata, mouseevent.ydata
        xmin, xmax, ymin, ymax = self.get_extent()
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        if ymin > ymax:
            ymin, ymax = ymax, ymin

        if x is not None and y is not None:
            inside = (xmin <= x <= xmax) and (ymin <= y <= ymax)
        else:
            inside = False

        return inside, {}

    def write_png(self, fname):
        """Write the image to png file with fname"""
        im = self.to_rgba(self._A[::-1] if self.origin == 'lower' else self._A,
                          bytes=True, norm=True)
        _png.write_png(im, fname)

    def set_data(self, A):
        """
        Set the image array.

        ACCEPTS: numpy/PIL Image A

        Note that this function does *not* update the normalization used.
        """
        # check if data is PIL Image without importing Image
        if hasattr(A, 'getpixel'):
            if A.mode == 'L':
                # greyscale image, but our logic assumes rgba:
                self._A = pil_to_array(A.convert('RGBA'))
            else:
                self._A = pil_to_array(A)
        else:
            self._A = cbook.safe_masked_invalid(A, copy=True)

        if (self._A.dtype != np.uint8 and
                not np.can_cast(self._A.dtype, float, "same_kind")):
            raise TypeError("Image data cannot be converted to float")

        if not (self._A.ndim == 2
                or self._A.ndim == 3 and self._A.shape[-1] in [3, 4]):
            raise TypeError("Invalid dimensions for image data")

        if self._A.ndim == 3:
            # If the input data has values outside the valid range (after
            # normalisation), we issue a warning and then clip X to the bounds
            # - otherwise casting wraps extreme values, hiding outliers and
            # making reliable interpretation impossible.
            high = 255 if np.issubdtype(self._A.dtype, np.integer) else 1
            if self._A.min() < 0 or high < self._A.max():
                _log.warning(
                    'Clipping input data to the valid range for imshow with '
                    'RGB data ([0..1] for floats or [0..255] for integers).'
                )
                self._A = np.clip(self._A, 0, high)
            # Cast unsupported integer types to uint8
            if self._A.dtype != np.uint8 and np.issubdtype(self._A.dtype,
                                                           np.integer):
                self._A = self._A.astype(np.uint8)

        self._imcache = None
        self._rgbacache = None
        self.stale = True

    def set_array(self, A):
        """
        Retained for backwards compatibility - use set_data instead

        ACCEPTS: numpy array A or PIL Image
        """
        # This also needs to be here to override the inherited
        # cm.ScalarMappable.set_array method so it is not invoked
        # by mistake.

        self.set_data(A)

    def get_interpolation(self):
        """
        Return the interpolation method the image uses when resizing.

        One of 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36',
        'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom',
        'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos', or 'none'.

        """
        return self._interpolation

    def set_interpolation(self, s):
        """
        Set the interpolation method the image uses when resizing.

        if None, use a value from rc setting. If 'none', the image is
        shown as is without interpolating. 'none' is only supported in
        agg, ps and pdf backends and will fall back to 'nearest' mode
        for other backends.

        .. ACCEPTS: ['nearest' | 'bilinear' | 'bicubic' | 'spline16' |
           'spline36' | 'hanning' | 'hamming' | 'hermite' | 'kaiser' |
           'quadric' | 'catrom' | 'gaussian' | 'bessel' | 'mitchell' |
           'sinc' | 'lanczos' | 'none' ]

        """
        if s is None:
            s = rcParams['image.interpolation']
        s = s.lower()
        if s not in _interpd_:
            raise ValueError('Illegal interpolation string')
        self._interpolation = s
        self.stale = True

    def can_composite(self):
        """
        Returns `True` if the image can be composited with its neighbors.
        """
        trans = self.get_transform()
        return (
            self._interpolation != 'none' and
            trans.is_affine and
            trans.is_separable)

    def set_resample(self, v):
        """
        Set whether or not image resampling is used.

        ACCEPTS: True|False
        """
        if v is None:
            v = rcParams['image.resample']
        self._resample = v
        self.stale = True

    def get_resample(self):
        """Return the image resample boolean."""
        return self._resample

    def set_filternorm(self, filternorm):
        """
        Set whether the resize filter norms the weights -- see
        help for imshow

        ACCEPTS: 0 or 1
        """
        if filternorm:
            self._filternorm = 1
        else:
            self._filternorm = 0

        self.stale = True

    def get_filternorm(self):
        """Return the filternorm setting."""
        return self._filternorm

    def set_filterrad(self, filterrad):
        """
        Set the resize filter radius only applicable to some
        interpolation schemes -- see help for imshow

        ACCEPTS: positive float
        """
        r = float(filterrad)
        if r <= 0:
            raise ValueError("The filter radius must be a positive number")
        self._filterrad = r
        self.stale = True

    def get_filterrad(self):
        """Return the filterrad setting."""
        return self._filterrad


class AxesImage(_ImageBase):
    def __str__(self):
        return "AxesImage(%g,%g;%gx%g)" % tuple(self.axes.bbox.bounds)

    def __init__(self, ax,
                 cmap=None,
                 norm=None,
                 interpolation=None,
                 origin=None,
                 extent=None,
                 filternorm=1,
                 filterrad=4.0,
                 resample=False,
                 **kwargs
                 ):

        """
        interpolation and cmap default to their rc settings

        cmap is a colors.Colormap instance
        norm is a colors.Normalize instance to map luminance to 0-1

        extent is data axes (left, right, bottom, top) for making image plots
        registered with data plots.  Default is to label the pixel
        centers with the zero-based row and column indices.

        Additional kwargs are matplotlib.artist properties

        """

        self._extent = extent

        super(AxesImage, self).__init__(
            ax,
            cmap=cmap,
            norm=norm,
            interpolation=interpolation,
            origin=origin,
            filternorm=filternorm,
            filterrad=filterrad,
            resample=resample,
            **kwargs
        )

    def get_window_extent(self, renderer=None):
        x0, x1, y0, y1 = self._extent
        bbox = Bbox.from_extents([x0, y0, x1, y1])
        return bbox.transformed(self.axes.transData)

    def make_image(self, renderer, magnification=1.0, unsampled=False):
        trans = self.get_transform()
        # image is created in the canvas coordinate.
        x1, x2, y1, y2 = self.get_extent()
        bbox = Bbox(np.array([[x1, y1], [x2, y2]]))
        transformed_bbox = TransformedBbox(bbox, trans)

        return self._make_image(
            self._A, bbox, transformed_bbox, self.axes.bbox, magnification,
            unsampled=unsampled)

    def _check_unsampled_image(self, renderer):
        """
        Return whether the image would be better drawn unsampled.
        """
        return (self.get_interpolation() == "none"
                and renderer.option_scale_image())

    def set_extent(self, extent):
        """
        extent is data axes (left, right, bottom, top) for making image plots

        This updates ax.dataLim, and, if autoscaling, sets viewLim
        to tightly fit the image, regardless of dataLim.  Autoscaling
        state is not changed, so following this with ax.autoscale_view
        will redo the autoscaling in accord with dataLim.
        """
        self._extent = xmin, xmax, ymin, ymax = extent
        corners = (xmin, ymin), (xmax, ymax)
        self.axes.update_datalim(corners)
        self.sticky_edges.x[:] = [xmin, xmax]
        self.sticky_edges.y[:] = [ymin, ymax]
        if self.axes._autoscaleXon:
            self.axes.set_xlim((xmin, xmax), auto=None)
        if self.axes._autoscaleYon:
            self.axes.set_ylim((ymin, ymax), auto=None)
        self.stale = True

    def get_extent(self):
        """Get the image extent: left, right, bottom, top"""
        if self._extent is not None:
            return self._extent
        else:
            sz = self.get_size()
            numrows, numcols = sz
            if self.origin == 'upper':
                return (-0.5, numcols-0.5, numrows-0.5, -0.5)
            else:
                return (-0.5, numcols-0.5, -0.5, numrows-0.5)

    def get_cursor_data(self, event):
        """Get the cursor data for a given event"""
        xmin, xmax, ymin, ymax = self.get_extent()
        if self.origin == 'upper':
            ymin, ymax = ymax, ymin
        arr = self.get_array()
        data_extent = Bbox([[ymin, xmin], [ymax, xmax]])
        array_extent = Bbox([[0, 0], arr.shape[:2]])
        trans = BboxTransform(boxin=data_extent, boxout=array_extent)
        y, x = event.ydata, event.xdata
        point = trans.transform_point([y, x])
        if any(np.isnan(point)):
            return None
        i, j = point.astype(int)
        # Clip the coordinates at array bounds
        if not (0 <= i < arr.shape[0]) or not (0 <= j < arr.shape[1]):
            return None
        else:
            return arr[i, j]


class NonUniformImage(AxesImage):
    def __init__(self, ax, **kwargs):
        """
        kwargs are identical to those for AxesImage, except
        that 'nearest' and 'bilinear' are the only supported 'interpolation'
        options.
        """
        interp = kwargs.pop('interpolation', 'nearest')
        super(NonUniformImage, self).__init__(ax, **kwargs)
        self.set_interpolation(interp)

    def _check_unsampled_image(self, renderer):
        """
        return False. Do not use unsampled image.
        """
        return False

    def make_image(self, renderer, magnification=1.0, unsampled=False):
        if self._A is None:
            raise RuntimeError('You must first set the image array')

        if unsampled:
            raise ValueError('unsampled not supported on NonUniformImage')

        A = self._A
        if A.ndim == 2:
            if A.dtype != np.uint8:
                A = self.to_rgba(A, bytes=True)
                self.is_grayscale = self.cmap.is_gray()
            else:
                A = np.repeat(A[:, :, np.newaxis], 4, 2)
                A[:, :, 3] = 255
                self.is_grayscale = True
        else:
            if A.dtype != np.uint8:
                A = (255*A).astype(np.uint8)
            if A.shape[2] == 3:
                B = np.zeros(tuple(list(A.shape[0:2]) + [4]), np.uint8)
                B[:, :, 0:3] = A
                B[:, :, 3] = 255
                A = B
            self.is_grayscale = False

        x0, y0, v_width, v_height = self.axes.viewLim.bounds
        l, b, r, t = self.axes.bbox.extents
        width = (np.round(r) + 0.5) - (np.round(l) - 0.5)
        height = (np.round(t) + 0.5) - (np.round(b) - 0.5)
        width *= magnification
        height *= magnification
        im = _image.pcolor(self._Ax, self._Ay, A,
                           int(height), int(width),
                           (x0, x0+v_width, y0, y0+v_height),
                           _interpd_[self._interpolation])

        return im, l, b, IdentityTransform()

    def set_data(self, x, y, A):
        """
        Set the grid for the pixel centers, and the pixel values.

          *x* and *y* are monotonic 1-D ndarrays of lengths N and M,
             respectively, specifying pixel centers

          *A* is an (M,N) ndarray or masked array of values to be
            colormapped, or a (M,N,3) RGB array, or a (M,N,4) RGBA
            array.
        """
        x = np.array(x, np.float32)
        y = np.array(y, np.float32)
        A = cbook.safe_masked_invalid(A, copy=True)
        if not (x.ndim == y.ndim == 1 and A.shape[0:2] == y.shape + x.shape):
            raise TypeError("Axes don't match array shape")
        if A.ndim not in [2, 3]:
            raise TypeError("Can only plot 2D or 3D data")
        if A.ndim == 3 and A.shape[2] not in [1, 3, 4]:
            raise TypeError("3D arrays must have three (RGB) "
                            "or four (RGBA) color components")
        if A.ndim == 3 and A.shape[2] == 1:
            A.shape = A.shape[0:2]
        self._A = A
        self._Ax = x
        self._Ay = y
        self._imcache = None

        self.stale = True

    def set_array(self, *args):
        raise NotImplementedError('Method not supported')

    def set_interpolation(self, s):
        """
        Parameters
        ----------
        s : str, None
            Either 'nearest', 'bilinear', or ``None``.
        """
        if s is not None and s not in ('nearest', 'bilinear'):
            raise NotImplementedError('Only nearest neighbor and '
                                      'bilinear interpolations are supported')
        AxesImage.set_interpolation(self, s)

    def get_extent(self):
        if self._A is None:
            raise RuntimeError('Must set data first')
        return self._Ax[0], self._Ax[-1], self._Ay[0], self._Ay[-1]

    def set_filternorm(self, s):
        pass

    def set_filterrad(self, s):
        pass

    def set_norm(self, norm):
        if self._A is not None:
            raise RuntimeError('Cannot change colors after loading data')
        super(NonUniformImage, self).set_norm(norm)

    def set_cmap(self, cmap):
        if self._A is not None:
            raise RuntimeError('Cannot change colors after loading data')
        super(NonUniformImage, self).set_cmap(cmap)


class PcolorImage(AxesImage):
    """
    Make a pcolor-style plot with an irregular rectangular grid.

    This uses a variation of the original irregular image code,
    and it is used by pcolorfast for the corresponding grid type.
    """
    def __init__(self, ax,
                 x=None,
                 y=None,
                 A=None,
                 cmap=None,
                 norm=None,
                 **kwargs
                 ):
        """
        cmap defaults to its rc setting

        cmap is a colors.Colormap instance
        norm is a colors.Normalize instance to map luminance to 0-1

        Additional kwargs are matplotlib.artist properties
        """
        super(PcolorImage, self).__init__(ax, norm=norm, cmap=cmap)
        self.update(kwargs)
        if A is not None:
            self.set_data(x, y, A)

    def make_image(self, renderer, magnification=1.0, unsampled=False):
        if self._A is None:
            raise RuntimeError('You must first set the image array')
        if unsampled:
            raise ValueError('unsampled not supported on PColorImage')
        fc = self.axes.patch.get_facecolor()
        bg = mcolors.to_rgba(fc, 0)
        bg = (np.array(bg)*255).astype(np.uint8)
        l, b, r, t = self.axes.bbox.extents
        width = (np.round(r) + 0.5) - (np.round(l) - 0.5)
        height = (np.round(t) + 0.5) - (np.round(b) - 0.5)
        # The extra cast-to-int is only needed for python2
        width = int(np.round(width * magnification))
        height = int(np.round(height * magnification))
        if self._rgbacache is None:
            A = self.to_rgba(self._A, bytes=True)
            self._rgbacache = A
            if self._A.ndim == 2:
                self.is_grayscale = self.cmap.is_gray()
        else:
            A = self._rgbacache
        vl = self.axes.viewLim
        im = _image.pcolor2(self._Ax, self._Ay, A,
                            height,
                            width,
                            (vl.x0, vl.x1, vl.y0, vl.y1),
                            bg)
        return im, l, b, IdentityTransform()

    def _check_unsampled_image(self, renderer):
        return False

    def set_data(self, x, y, A):
        """
        Set the grid for the rectangle boundaries, and the data values.

          *x* and *y* are monotonic 1-D ndarrays of lengths N+1 and M+1,
             respectively, specifying rectangle boundaries.  If None,
             they will be created as uniform arrays from 0 through N
             and 0 through M, respectively.

          *A* is an (M,N) ndarray or masked array of values to be
            colormapped, or a (M,N,3) RGB array, or a (M,N,4) RGBA
            array.

        """
        A = cbook.safe_masked_invalid(A, copy=True)
        if x is None:
            x = np.arange(0, A.shape[1]+1, dtype=np.float64)
        else:
            x = np.array(x, np.float64).ravel()
        if y is None:
            y = np.arange(0, A.shape[0]+1, dtype=np.float64)
        else:
            y = np.array(y, np.float64).ravel()

        if A.shape[:2] != (y.size-1, x.size-1):
            raise ValueError(
                "Axes don't match array shape. Got %s, expected %s." %
                (A.shape[:2], (y.size - 1, x.size - 1)))
        if A.ndim not in [2, 3]:
            raise ValueError("A must be 2D or 3D")
        if A.ndim == 3 and A.shape[2] == 1:
            A.shape = A.shape[:2]
        self.is_grayscale = False
        if A.ndim == 3:
            if A.shape[2] in [3, 4]:
                if ((A[:, :, 0] == A[:, :, 1]).all() and
                        (A[:, :, 0] == A[:, :, 2]).all()):
                    self.is_grayscale = True
            else:
                raise ValueError("3D arrays must have RGB or RGBA as last dim")

        # For efficient cursor readout, ensure x and y are increasing.
        if x[-1] < x[0]:
            x = x[::-1]
            A = A[:, ::-1]
        if y[-1] < y[0]:
            y = y[::-1]
            A = A[::-1]

        self._A = A
        self._Ax = x
        self._Ay = y
        self._rgbacache = None
        self.stale = True

    def set_array(self, *args):
        raise NotImplementedError('Method not supported')

    def get_cursor_data(self, event):
        """Get the cursor data for a given event"""
        x, y = event.xdata, event.ydata
        if (x < self._Ax[0] or x > self._Ax[-1] or
                y < self._Ay[0] or y > self._Ay[-1]):
            return None
        j = np.searchsorted(self._Ax, x) - 1
        i = np.searchsorted(self._Ay, y) - 1
        try:
            return self._A[i, j]
        except IndexError:
            return None


class FigureImage(_ImageBase):
    zorder = 0

    _interpolation = 'nearest'

    def __init__(self, fig,
                 cmap=None,
                 norm=None,
                 offsetx=0,
                 offsety=0,
                 origin=None,
                 **kwargs
                 ):
        """
        cmap is a colors.Colormap instance
        norm is a colors.Normalize instance to map luminance to 0-1

        kwargs are an optional list of Artist keyword args
        """
        super(FigureImage, self).__init__(
            None,
            norm=norm,
            cmap=cmap,
            origin=origin
        )
        self.figure = fig
        self.ox = offsetx
        self.oy = offsety
        self.update(kwargs)
        self.magnification = 1.0

    def get_extent(self):
        """Get the image extent: left, right, bottom, top"""
        numrows, numcols = self.get_size()
        return (-0.5 + self.ox, numcols-0.5 + self.ox,
                -0.5 + self.oy, numrows-0.5 + self.oy)

    def make_image(self, renderer, magnification=1.0, unsampled=False):
        fac = renderer.dpi/self.figure.dpi
        # fac here is to account for pdf, eps, svg backends where
        # figure.dpi is set to 72.  This means we need to scale the
        # image (using magification) and offset it appropriately.
        bbox = Bbox([[self.ox/fac, self.oy/fac],
                     [(self.ox/fac + self._A.shape[1]),
                     (self.oy/fac + self._A.shape[0])]])
        width, height = self.figure.get_size_inches()
        width *= renderer.dpi
        height *= renderer.dpi
        clip = Bbox([[0, 0], [width, height]])

        return self._make_image(
            self._A, bbox, bbox, clip, magnification=magnification / fac,
            unsampled=unsampled, round_to_pixel_border=False)

    def set_data(self, A):
        """Set the image array."""
        cm.ScalarMappable.set_array(self,
                                    cbook.safe_masked_invalid(A, copy=True))
        self.stale = True


class BboxImage(_ImageBase):
    """The Image class whose size is determined by the given bbox."""
    def __init__(self, bbox,
                 cmap=None,
                 norm=None,
                 interpolation=None,
                 origin=None,
                 filternorm=1,
                 filterrad=4.0,
                 resample=False,
                 interp_at_native=True,
                 **kwargs
                 ):
        """
        cmap is a colors.Colormap instance
        norm is a colors.Normalize instance to map luminance to 0-1

        interp_at_native is a flag that determines whether or not
        interpolation should still be applied when the image is
        displayed at its native resolution.  A common use case for this
        is when displaying an image for annotational purposes; it is
        treated similarly to Photoshop (interpolation is only used when
        displaying the image at non-native resolutions).


        kwargs are an optional list of Artist keyword args

        """
        super(BboxImage, self).__init__(
            None,
            cmap=cmap,
            norm=norm,
            interpolation=interpolation,
            origin=origin,
            filternorm=filternorm,
            filterrad=filterrad,
            resample=resample,
            **kwargs
        )

        self.bbox = bbox
        self.interp_at_native = interp_at_native
        self._transform = IdentityTransform()

    def get_transform(self):
        return self._transform

    def get_window_extent(self, renderer=None):
        if renderer is None:
            renderer = self.get_figure()._cachedRenderer

        if isinstance(self.bbox, BboxBase):
            return self.bbox
        elif callable(self.bbox):
            return self.bbox(renderer)
        else:
            raise ValueError("unknown type of bbox")

    def contains(self, mouseevent):
        """Test whether the mouse event occurred within the image."""
        if callable(self._contains):
            return self._contains(self, mouseevent)

        if not self.get_visible():  # or self.get_figure()._renderer is None:
            return False, {}

        x, y = mouseevent.x, mouseevent.y
        inside = self.get_window_extent().contains(x, y)

        return inside, {}

    def make_image(self, renderer, magnification=1.0, unsampled=False):
        width, height = renderer.get_canvas_width_height()

        bbox_in = self.get_window_extent(renderer).frozen()
        bbox_in._points /= [width, height]
        bbox_out = self.get_window_extent(renderer)
        clip = Bbox([[0, 0], [width, height]])
        self._transform = BboxTransform(Bbox([[0, 0], [1, 1]]), clip)

        return self._make_image(
            self._A,
            bbox_in, bbox_out, clip, magnification, unsampled=unsampled)


def imread(fname, format=None):
    """
    Read an image from a file into an array.

    *fname* may be a string path, a valid URL, or a Python
    file-like object.  If using a file object, it must be opened in binary
    mode.

    If *format* is provided, will try to read file of that type,
    otherwise the format is deduced from the filename.  If nothing can
    be deduced, PNG is tried.

    Return value is a :class:`numpy.array`.  For grayscale images, the
    return array is MxN.  For RGB images, the return value is MxNx3.
    For RGBA images the return value is MxNx4.

    matplotlib can only read PNGs natively, but if `PIL
    <http://www.pythonware.com/products/pil/>`_ is installed, it will
    use it to load the image and return an array (if possible) which
    can be used with :func:`~matplotlib.pyplot.imshow`. Note, URL strings
    may not be compatible with PIL. Check the PIL documentation for more
    information.
    """

    def pilread(fname):
        """try to load the image with PIL or return None"""
        try:
            from PIL import Image
        except ImportError:
            return None
        with Image.open(fname) as image:
            return pil_to_array(image)

    handlers = {'png': _png.read_png, }
    if format is None:
        if isinstance(fname, six.string_types):
            parsed = urlparse(fname)
            # If the string is a URL, assume png
            if len(parsed.scheme) > 1:
                ext = 'png'
            else:
                basename, ext = os.path.splitext(fname)
                ext = ext.lower()[1:]
        elif hasattr(fname, 'name'):
            basename, ext = os.path.splitext(fname.name)
            ext = ext.lower()[1:]
        else:
            ext = 'png'
    else:
        ext = format

    if ext not in handlers:
        im = pilread(fname)
        if im is None:
            raise ValueError('Only know how to handle extensions: %s; '
                             'with Pillow installed matplotlib can handle '
                             'more images' % list(handlers))
        return im

    handler = handlers[ext]

    # To handle Unicode filenames, we pass a file object to the PNG
    # reader extension, since Python handles them quite well, but it's
    # tricky in C.
    if isinstance(fname, six.string_types):
        parsed = urlparse(fname)
        # If fname is a URL, download the data
        if len(parsed.scheme) > 1:
            fd = BytesIO(urlopen(fname).read())
            return handler(fd)
        else:
            with open(fname, 'rb') as fd:
                return handler(fd)
    else:
        return handler(fname)


def imsave(fname, arr, vmin=None, vmax=None, cmap=None, format=None,
           origin=None, dpi=100):
    """
    Save an array as in image file.

    The output formats available depend on the backend being used.

    Parameters
    ----------
    fname : str or file-like
        Path string to a filename, or a Python file-like object.
        If *format* is *None* and *fname* is a string, the output
        format is deduced from the extension of the filename.
    arr : array-like
        An MxN (luminance), MxNx3 (RGB) or MxNx4 (RGBA) array.
    vmin, vmax: [ None | scalar ]
        *vmin* and *vmax* set the color scaling for the image by fixing the
        values that map to the colormap color limits. If either *vmin*
        or *vmax* is None, that limit is determined from the *arr*
        min/max value.
    cmap : matplotlib.colors.Colormap, optional
        For example, ``cm.viridis``.  If ``None``, defaults to the
        ``image.cmap`` rcParam.
    format : str
        One of the file extensions supported by the active backend.  Most
        backends support png, pdf, ps, eps and svg.
    origin : [ 'upper' | 'lower' ]
        Indicates whether the ``(0, 0)`` index of the array is in the
        upper left or lower left corner of the axes.  Defaults to the
        ``image.origin`` rcParam.
    dpi : int
        The DPI to store in the metadata of the file.  This does not affect the
        resolution of the output image.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    if isinstance(fname, getattr(os, "PathLike", ())):
        fname = os.fspath(fname)
    if (format == 'png'
        or (format is None
            and isinstance(fname, six.string_types)
            and fname.lower().endswith('.png'))):
        image = AxesImage(None, cmap=cmap, origin=origin)
        image.set_data(arr)
        image.set_clim(vmin, vmax)
        image.write_png(fname)
    else:
        fig = Figure(dpi=dpi, frameon=False)
        FigureCanvas(fig)
        fig.figimage(arr, cmap=cmap, vmin=vmin, vmax=vmax, origin=origin,
                     resize=True)
        fig.savefig(fname, dpi=dpi, format=format, transparent=True)


def pil_to_array(pilImage):
    """Load a PIL image and return it as a numpy array.

    Grayscale images are returned as ``(M, N)`` arrays.  RGB images are
    returned as ``(M, N, 3)`` arrays.  RGBA images are returned as ``(M, N,
    4)`` arrays.
    """
    if pilImage.mode in ['RGBA', 'RGBX', 'RGB', 'L']:
        # return MxNx4 RGBA, MxNx3 RBA, or MxN luminance array
        return np.asarray(pilImage)
    elif pilImage.mode.startswith('I;16'):
        # return MxN luminance array of uint16
        raw = pilImage.tobytes('raw', pilImage.mode)
        if pilImage.mode.endswith('B'):
            x = np.fromstring(raw, '>u2')
        else:
            x = np.fromstring(raw, '<u2')
        return x.reshape(pilImage.size[::-1]).astype('=u2')
    else:  # try to convert to an rgba image
        try:
            pilImage = pilImage.convert('RGBA')
        except ValueError:
            raise RuntimeError('Unknown image mode')
        return np.asarray(pilImage)  # return MxNx4 RGBA array


def thumbnail(infile, thumbfile, scale=0.1, interpolation='bilinear',
              preview=False):
    """
    make a thumbnail of image in *infile* with output filename
    *thumbfile*.

      *infile* the image file -- must be PNG or Pillow-readable if you
         have `Pillow <http://python-pillow.org/>`_ installed

      *thumbfile*
        the thumbnail filename

      *scale*
        the scale factor for the thumbnail

      *interpolation*
        the interpolation scheme used in the resampling


      *preview*
        if True, the default backend (presumably a user interface
        backend) will be used which will cause a figure to be raised
        if :func:`~matplotlib.pyplot.show` is called.  If it is False,
        a pure image backend will be used depending on the extension,
        'png'->FigureCanvasAgg, 'pdf'->FigureCanvasPdf,
        'svg'->FigureCanvasSVG


    See examples/misc/image_thumbnail.py.

    .. htmlonly::

        :ref:`sphx_glr_gallery_misc_image_thumbnail_sgskip.py`

    Return value is the figure instance containing the thumbnail

    """
    basedir, basename = os.path.split(infile)
    baseout, extout = os.path.splitext(thumbfile)

    im = imread(infile)
    rows, cols, depth = im.shape

    # this doesn't really matter, it will cancel in the end, but we
    # need it for the mpl API
    dpi = 100

    height = rows / dpi * scale
    width = cols / dpi * scale

    extension = extout.lower()

    if preview:
        # let the UI backend do everything
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(width, height), dpi=dpi)
    else:
        if extension == '.png':
            from matplotlib.backends.backend_agg \
                import FigureCanvasAgg as FigureCanvas
        elif extension == '.pdf':
            from matplotlib.backends.backend_pdf \
                import FigureCanvasPdf as FigureCanvas
        elif extension == '.svg':
            from matplotlib.backends.backend_svg \
                import FigureCanvasSVG as FigureCanvas
        else:
            raise ValueError("Can only handle "
                             "extensions 'png', 'svg' or 'pdf'")

        from matplotlib.figure import Figure
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas(fig)

    ax = fig.add_axes([0, 0, 1, 1], aspect='auto',
                      frameon=False, xticks=[], yticks=[])

    basename, ext = os.path.splitext(basename)
    ax.imshow(im, aspect='auto', resample=True, interpolation=interpolation)
    fig.savefig(thumbfile, dpi=dpi)
    return fig
