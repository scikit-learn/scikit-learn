from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
from .axes_divider import make_axes_locatable, Size, locatable_axes_factory
import sys
from .mpl_axes import Axes


def make_rgb_axes(ax, pad=0.01, axes_class=None, add_all=True):
    """
    pad : fraction of the axes height.
    """

    divider = make_axes_locatable(ax)

    pad_size = Size.Fraction(pad, Size.AxesY(ax))

    xsize = Size.Fraction((1.-2.*pad)/3., Size.AxesX(ax))
    ysize = Size.Fraction((1.-2.*pad)/3., Size.AxesY(ax))

    divider.set_horizontal([Size.AxesX(ax), pad_size, xsize])
    divider.set_vertical([ysize, pad_size, ysize, pad_size, ysize])

    ax.set_axes_locator(divider.new_locator(0, 0, ny1=-1))

    ax_rgb = []
    if axes_class is None:
        try:
            axes_class = locatable_axes_factory(ax._axes_class)
        except AttributeError:
            axes_class = locatable_axes_factory(type(ax))

    for ny in [4, 2, 0]:
        ax1 = axes_class(ax.get_figure(),
                         ax.get_position(original=True),
                         sharex=ax, sharey=ax)
        locator = divider.new_locator(nx=2, ny=ny)
        ax1.set_axes_locator(locator)
        for t in ax1.yaxis.get_ticklabels() + ax1.xaxis.get_ticklabels():
            t.set_visible(False)
        try:
            for axis in ax1.axis.values():
                axis.major_ticklabels.set_visible(False)
        except AttributeError:
            pass

        ax_rgb.append(ax1)

    if add_all:
        fig = ax.get_figure()
        for ax1 in ax_rgb:
            fig.add_axes(ax1)

    return ax_rgb


def imshow_rgb(ax, r, g, b, **kwargs):
    ny, nx = r.shape
    R = np.zeros([ny, nx, 3], dtype="d")
    R[:,:,0] = r
    G = np.zeros_like(R)
    G[:,:,1] = g
    B = np.zeros_like(R)
    B[:,:,2] = b

    RGB = R + G + B

    im_rgb = ax.imshow(RGB, **kwargs)

    return im_rgb


class RGBAxesBase(object):
    """base class for a 4-panel imshow (RGB, R, G, B)

    Layout:
    +---------------+-----+
    |               |  R  |
    +               +-----+
    |      RGB      |  G  |
    +               +-----+
    |               |  B  |
    +---------------+-----+

    Attributes
    ----------
    _defaultAxesClass : matplotlib.axes.Axes
        defaults to 'Axes' in RGBAxes child class.
        No default in abstract base class
    RGB : _defaultAxesClass
        The axes object for the three-channel imshow
    R : _defaultAxesClass
        The axes object for the red channel imshow
    G : _defaultAxesClass
        The axes object for the green channel imshow
    B : _defaultAxesClass
        The axes object for the blue channel imshow
    """
    def __init__(self, *kl, **kwargs):
        """
        Parameters
        ----------
        pad : float
            fraction of the axes height to put as padding.
            defaults to 0.0
        add_all : bool
            True: Add the {rgb, r, g, b} axes to the figure
            defaults to True.
        axes_class : matplotlib.axes.Axes

        kl :
            Unpacked into axes_class() init for RGB
        kwargs :
            Unpacked into axes_class() init for RGB, R, G, B axes
        """
        pad = kwargs.pop("pad", 0.0)
        add_all = kwargs.pop("add_all", True)
        try:
            axes_class = kwargs.pop("axes_class", self._defaultAxesClass)
        except AttributeError:
            new_msg = ("A subclass of RGBAxesBase must have a "
                       "_defaultAxesClass attribute. If you are not sure which "
                       "axes class to use, consider using "
                       "mpl_toolkits.axes_grid1.mpl_axes.Axes.")
            six.reraise(AttributeError, AttributeError(new_msg),
                        sys.exc_info()[2])

        ax = axes_class(*kl, **kwargs)

        divider = make_axes_locatable(ax)

        pad_size = Size.Fraction(pad, Size.AxesY(ax))

        xsize = Size.Fraction((1.-2.*pad)/3., Size.AxesX(ax))
        ysize = Size.Fraction((1.-2.*pad)/3., Size.AxesY(ax))

        divider.set_horizontal([Size.AxesX(ax), pad_size, xsize])
        divider.set_vertical([ysize, pad_size, ysize, pad_size, ysize])

        ax.set_axes_locator(divider.new_locator(0, 0, ny1=-1))

        ax_rgb = []
        for ny in [4, 2, 0]:
            ax1 = axes_class(ax.get_figure(),
                             ax.get_position(original=True),
                             sharex=ax, sharey=ax, **kwargs)
            locator = divider.new_locator(nx=2, ny=ny)
            ax1.set_axes_locator(locator)
            ax1.axis[:].toggle(ticklabels=False)
            ax_rgb.append(ax1)

        self.RGB = ax
        self.R, self.G, self.B = ax_rgb

        if add_all:
            fig = ax.get_figure()
            fig.add_axes(ax)
            self.add_RGB_to_figure()

        self._config_axes()

    def _config_axes(self, line_color='w', marker_edge_color='w'):
        """Set the line color and ticks for the axes

        Parameters
        ----------
        line_color : any matplotlib color
        marker_edge_color : any matplotlib color
        """
        for ax1 in [self.RGB, self.R, self.G, self.B]:
            ax1.axis[:].line.set_color(line_color)
            ax1.axis[:].major_ticks.set_markeredgecolor(marker_edge_color)

    def add_RGB_to_figure(self):
        """Add the red, green and blue axes to the RGB composite's axes figure
        """
        self.RGB.get_figure().add_axes(self.R)
        self.RGB.get_figure().add_axes(self.G)
        self.RGB.get_figure().add_axes(self.B)

    def imshow_rgb(self, r, g, b, **kwargs):
        """Create the four images {rgb, r, g, b}

        Parameters
        ----------
        r : array-like
            The red array
        g : array-like
            The green array
        b : array-like
            The blue array
        kwargs : imshow kwargs
            kwargs get unpacked into the imshow calls for the four images

        Returns
        -------
        rgb : matplotlib.image.AxesImage
        r : matplotlib.image.AxesImage
        g : matplotlib.image.AxesImage
        b : matplotlib.image.AxesImage
        """
        if not (r.shape == g.shape == b.shape):
            raise ValueError('Input shapes do not match.'
                             '\nr.shape = {}'
                             '\ng.shape = {}'
                             '\nb.shape = {}'
                             .format(r.shape, g.shape, b.shape))
        RGB = np.dstack([r, g, b])
        R = np.zeros_like(RGB)
        R[:,:,0] = r
        G = np.zeros_like(RGB)
        G[:,:,1] = g
        B = np.zeros_like(RGB)
        B[:,:,2] = b

        im_rgb = self.RGB.imshow(RGB, **kwargs)
        im_r = self.R.imshow(R, **kwargs)
        im_g = self.G.imshow(G, **kwargs)
        im_b = self.B.imshow(B, **kwargs)

        return im_rgb, im_r, im_g, im_b


class RGBAxes(RGBAxesBase):
    _defaultAxesClass = Axes
