"""
Builtin colormaps, colormap handling utilities, and the `ScalarMappable` mixin.

See :doc:`/gallery/color/colormap_reference` for a list of builtin colormaps.
See :doc:`/tutorials/colors/colormaps` for an in-depth discussion of colormaps.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
from numpy import ma
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib._cm import datad
from matplotlib._cm_listed import cmaps as cmaps_listed


cmap_d = {}


# reverse all the colormaps.
# reversed colormaps have '_r' appended to the name.


def _reverser(f):
    def freversed(x):
        return f(1 - x)
    return freversed


def revcmap(data):
    """Can only handle specification *data* in dictionary format."""
    data_r = {}
    for key, val in six.iteritems(data):
        if callable(val):
            valnew = _reverser(val)
            # This doesn't work: lambda x: val(1-x)
            # The same "val" (the first one) is used
            # each time, so the colors are identical
            # and the result is shades of gray.
        else:
            # Flip x and exchange the y values facing x = 0 and x = 1.
            valnew = [(1.0 - x, y1, y0) for x, y0, y1 in reversed(val)]
        data_r[key] = valnew
    return data_r


def _reverse_cmap_spec(spec):
    """Reverses cmap specification *spec*, can handle both dict and tuple
    type specs."""

    if 'listed' in spec:
        return {'listed': spec['listed'][::-1]}

    if 'red' in spec:
        return revcmap(spec)
    else:
        revspec = list(reversed(spec))
        if len(revspec[0]) == 2:    # e.g., (1, (1.0, 0.0, 1.0))
            revspec = [(1.0 - a, b) for a, b in revspec]
        return revspec


def _generate_cmap(name, lutsize):
    """Generates the requested cmap from its *name*.  The lut size is
    *lutsize*."""

    spec = datad[name]

    # Generate the colormap object.
    if 'red' in spec:
        return colors.LinearSegmentedColormap(name, spec, lutsize)
    elif 'listed' in spec:
        return colors.ListedColormap(spec['listed'], name)
    else:
        return colors.LinearSegmentedColormap.from_list(name, spec, lutsize)

LUTSIZE = mpl.rcParams['image.lut']

# Generate the reversed specifications (all at once, to avoid
# modify-when-iterating).
datad.update({cmapname + '_r': _reverse_cmap_spec(spec)
              for cmapname, spec in six.iteritems(datad)})

# Precache the cmaps with ``lutsize = LUTSIZE``.
# Also add the reversed ones added in the section above:
for cmapname in datad:
    cmap_d[cmapname] = _generate_cmap(cmapname, LUTSIZE)

cmap_d.update(cmaps_listed)

locals().update(cmap_d)


# Continue with definitions ...


def register_cmap(name=None, cmap=None, data=None, lut=None):
    """
    Add a colormap to the set recognized by :func:`get_cmap`.

    It can be used in two ways::

        register_cmap(name='swirly', cmap=swirly_cmap)

        register_cmap(name='choppy', data=choppydata, lut=128)

    In the first case, *cmap* must be a :class:`matplotlib.colors.Colormap`
    instance.  The *name* is optional; if absent, the name will
    be the :attr:`~matplotlib.colors.Colormap.name` attribute of the *cmap*.

    In the second case, the three arguments are passed to
    the :class:`~matplotlib.colors.LinearSegmentedColormap` initializer,
    and the resulting colormap is registered.

    """
    if name is None:
        try:
            name = cmap.name
        except AttributeError:
            raise ValueError("Arguments must include a name or a Colormap")

    if not isinstance(name, six.string_types):
        raise ValueError("Colormap name must be a string")

    if isinstance(cmap, colors.Colormap):
        cmap_d[name] = cmap
        return

    # For the remainder, let exceptions propagate.
    if lut is None:
        lut = mpl.rcParams['image.lut']
    cmap = colors.LinearSegmentedColormap(name, data, lut)
    cmap_d[name] = cmap


def get_cmap(name=None, lut=None):
    """
    Get a colormap instance, defaulting to rc values if *name* is None.

    Colormaps added with :func:`register_cmap` take precedence over
    built-in colormaps.

    If *name* is a :class:`matplotlib.colors.Colormap` instance, it will be
    returned.

    If *lut* is not None it must be an integer giving the number of
    entries desired in the lookup table, and *name* must be a standard
    mpl colormap name.
    """
    if name is None:
        name = mpl.rcParams['image.cmap']

    if isinstance(name, colors.Colormap):
        return name

    if name in cmap_d:
        if lut is None:
            return cmap_d[name]
        else:
            return cmap_d[name]._resample(lut)
    else:
        raise ValueError(
            "Colormap %s is not recognized. Possible values are: %s"
            % (name, ', '.join(sorted(cmap_d))))


class ScalarMappable(object):
    """
    This is a mixin class to support scalar data to RGBA mapping.
    The ScalarMappable makes use of data normalization before returning
    RGBA colors from the given colormap.

    """
    def __init__(self, norm=None, cmap=None):
        r"""

        Parameters
        ----------
        norm : :class:`matplotlib.colors.Normalize` instance
            The normalizing object which scales data, typically into the
            interval ``[0, 1]``.
            If *None*, *norm* defaults to a *colors.Normalize* object which
            initializes its scaling based on the first data processed.
        cmap : str or :class:`~matplotlib.colors.Colormap` instance
            The colormap used to map normalized data values to RGBA colors.
        """

        self.callbacksSM = cbook.CallbackRegistry()

        if cmap is None:
            cmap = get_cmap()
        if norm is None:
            norm = colors.Normalize()

        self._A = None
        #: The Normalization instance of this ScalarMappable.
        self.norm = norm
        #: The Colormap instance of this ScalarMappable.
        self.cmap = get_cmap(cmap)
        #: The last colorbar associated with this ScalarMappable. May be None.
        self.colorbar = None
        self.update_dict = {'array': False}

    def to_rgba(self, x, alpha=None, bytes=False, norm=True):
        """
        Return a normalized rgba array corresponding to *x*.

        In the normal case, *x* is a 1-D or 2-D sequence of scalars, and
        the corresponding ndarray of rgba values will be returned,
        based on the norm and colormap set for this ScalarMappable.

        There is one special case, for handling images that are already
        rgb or rgba, such as might have been read from an image file.
        If *x* is an ndarray with 3 dimensions,
        and the last dimension is either 3 or 4, then it will be
        treated as an rgb or rgba array, and no mapping will be done.
        The array can be uint8, or it can be floating point with
        values in the 0-1 range; otherwise a ValueError will be raised.
        If it is a masked array, the mask will be ignored.
        If the last dimension is 3, the *alpha* kwarg (defaulting to 1)
        will be used to fill in the transparency.  If the last dimension
        is 4, the *alpha* kwarg is ignored; it does not
        replace the pre-existing alpha.  A ValueError will be raised
        if the third dimension is other than 3 or 4.

        In either case, if *bytes* is *False* (default), the rgba
        array will be floats in the 0-1 range; if it is *True*,
        the returned rgba array will be uint8 in the 0 to 255 range.

        If norm is False, no normalization of the input data is
        performed, and it is assumed to be in the range (0-1).

        """
        # First check for special case, image input:
        try:
            if x.ndim == 3:
                if x.shape[2] == 3:
                    if alpha is None:
                        alpha = 1
                    if x.dtype == np.uint8:
                        alpha = np.uint8(alpha * 255)
                    m, n = x.shape[:2]
                    xx = np.empty(shape=(m, n, 4), dtype=x.dtype)
                    xx[:, :, :3] = x
                    xx[:, :, 3] = alpha
                elif x.shape[2] == 4:
                    xx = x
                else:
                    raise ValueError("third dimension must be 3 or 4")
                if xx.dtype.kind == 'f':
                    if norm and xx.max() > 1 or xx.min() < 0:
                        raise ValueError("Floating point image RGB values "
                                         "must be in the 0..1 range.")
                    if bytes:
                        xx = (xx * 255).astype(np.uint8)
                elif xx.dtype == np.uint8:
                    if not bytes:
                        xx = xx.astype(np.float32) / 255
                else:
                    raise ValueError("Image RGB array must be uint8 or "
                                     "floating point; found %s" % xx.dtype)
                return xx
        except AttributeError:
            # e.g., x is not an ndarray; so try mapping it
            pass

        # This is the normal case, mapping a scalar array:
        x = ma.asarray(x)
        if norm:
            x = self.norm(x)
        rgba = self.cmap(x, alpha=alpha, bytes=bytes)
        return rgba

    def set_array(self, A):
        """Set the image array from numpy array *A*.

        .. ACCEPTS: ndarray

        Parameters
        ----------
        A : ndarray
        """
        self._A = A
        self.update_dict['array'] = True

    def get_array(self):
        'Return the array'
        return self._A

    def get_cmap(self):
        'return the colormap'
        return self.cmap

    def get_clim(self):
        'return the min, max of the color limits for image scaling'
        return self.norm.vmin, self.norm.vmax

    def set_clim(self, vmin=None, vmax=None):
        """
        set the norm limits for image scaling; if *vmin* is a length2
        sequence, interpret it as ``(vmin, vmax)`` which is used to
        support setp

        ACCEPTS: a length 2 sequence of floats; may be overridden in methods
        that have ``vmin`` and ``vmax`` kwargs.
        """
        if vmax is None:
            try:
                vmin, vmax = vmin
            except (TypeError, ValueError):
                pass
        if vmin is not None:
            self.norm.vmin = colors._sanitize_extrema(vmin)
        if vmax is not None:
            self.norm.vmax = colors._sanitize_extrema(vmax)
        self.changed()

    def set_cmap(self, cmap):
        """
        set the colormap for luminance data

        ACCEPTS: a colormap or registered colormap name
        """
        cmap = get_cmap(cmap)
        self.cmap = cmap
        self.changed()

    def set_norm(self, norm):
        """Set the normalization instance.

        .. ACCEPTS: `.Normalize`

        Parameters
        ----------
        norm : `.Normalize`
        """
        if norm is None:
            norm = colors.Normalize()
        self.norm = norm
        self.changed()

    def autoscale(self):
        """
        Autoscale the scalar limits on the norm instance using the
        current array
        """
        if self._A is None:
            raise TypeError('You must first set_array for mappable')
        self.norm.autoscale(self._A)
        self.changed()

    def autoscale_None(self):
        """
        Autoscale the scalar limits on the norm instance using the
        current array, changing only limits that are None
        """
        if self._A is None:
            raise TypeError('You must first set_array for mappable')
        self.norm.autoscale_None(self._A)
        self.changed()

    def add_checker(self, checker):
        """
        Add an entry to a dictionary of boolean flags
        that are set to True when the mappable is changed.
        """
        self.update_dict[checker] = False

    def check_update(self, checker):
        """
        If mappable has changed since the last check,
        return True; else return False
        """
        if self.update_dict[checker]:
            self.update_dict[checker] = False
            return True
        return False

    def changed(self):
        """
        Call this whenever the mappable is changed to notify all the
        callbackSM listeners to the 'changed' signal
        """
        self.callbacksSM.process('changed', self)

        for key in self.update_dict:
            self.update_dict[key] = True
        self.stale = True
