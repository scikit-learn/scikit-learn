"""
This module is to support *bbox_inches* option in savefig command.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from matplotlib.transforms import Bbox, TransformedBbox, Affine2D


def adjust_bbox(fig, bbox_inches, fixed_dpi=None):
    """
    Temporarily adjust the figure so that only the specified area
    (bbox_inches) is saved.

    It modifies fig.bbox, fig.bbox_inches,
    fig.transFigure._boxout, and fig.patch.  While the figure size
    changes, the scale of the original figure is conserved.  A
    function which restores the original values are returned.
    """

    origBbox = fig.bbox
    origBboxInches = fig.bbox_inches
    _boxout = fig.transFigure._boxout

    asp_list = []
    locator_list = []
    for ax in fig.axes:
        pos = ax.get_position(original=False).frozen()
        locator_list.append(ax.get_axes_locator())
        asp_list.append(ax.get_aspect())

        def _l(a, r, pos=pos):
            return pos
        ax.set_axes_locator(_l)
        ax.set_aspect("auto")

    def restore_bbox():

        for ax, asp, loc in zip(fig.axes, asp_list, locator_list):
            ax.set_aspect(asp)
            ax.set_axes_locator(loc)

        fig.bbox = origBbox
        fig.bbox_inches = origBboxInches
        fig.transFigure._boxout = _boxout
        fig.transFigure.invalidate()
        fig.patch.set_bounds(0, 0, 1, 1)

    if fixed_dpi is not None:
        tr = Affine2D().scale(fixed_dpi)
        dpi_scale = fixed_dpi / fig.dpi
    else:
        tr = Affine2D().scale(fig.dpi)
        dpi_scale = 1.

    _bbox = TransformedBbox(bbox_inches, tr)

    fig.bbox_inches = Bbox.from_bounds(0, 0,
                                       bbox_inches.width, bbox_inches.height)
    x0, y0 = _bbox.x0, _bbox.y0
    w1, h1 = fig.bbox.width * dpi_scale, fig.bbox.height * dpi_scale
    fig.transFigure._boxout = Bbox.from_bounds(-x0, -y0, w1, h1)
    fig.transFigure.invalidate()

    fig.bbox = TransformedBbox(fig.bbox_inches, tr)

    fig.patch.set_bounds(x0 / w1, y0 / h1,
                         fig.bbox.width / w1, fig.bbox.height / h1)

    return restore_bbox


def process_figure_for_rasterizing(fig, bbox_inches_restore, fixed_dpi=None):
    """
    This need to be called when figure dpi changes during the drawing
    (e.g., rasterizing). It recovers the bbox and re-adjust it with
    the new dpi.
    """

    bbox_inches, restore_bbox = bbox_inches_restore
    restore_bbox()
    r = adjust_bbox(fig, bbox_inches, fixed_dpi)

    return bbox_inches, r
