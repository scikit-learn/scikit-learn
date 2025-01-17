"""
Utility Routines for Working with Matplotlib Objects
====================================================
"""
import itertools
import io
import base64

import numpy as np

import warnings

import matplotlib
from matplotlib.colors import colorConverter
from matplotlib.path import Path
from matplotlib.markers import MarkerStyle
from matplotlib.transforms import Affine2D
from matplotlib import ticker


def export_color(color):
    """Convert matplotlib color code to hex color or RGBA color"""
    if color is None or colorConverter.to_rgba(color)[3] == 0:
        return "none"
    elif colorConverter.to_rgba(color)[3] == 1:
        rgb = colorConverter.to_rgb(color)
        return "#{0:02X}{1:02X}{2:02X}".format(*(int(255 * c) for c in rgb))
    else:
        c = colorConverter.to_rgba(color)
        return (
            "rgba("
            + ", ".join(str(int(np.round(val * 255))) for val in c[:3])
            + ", "
            + str(c[3])
            + ")"
        )


def _many_to_one(input_dict):
    """Convert a many-to-one mapping to a one-to-one mapping"""
    return dict((key, val) for keys, val in input_dict.items() for key in keys)


LINESTYLES = _many_to_one(
    {
        ("solid", "-", (None, None)): "none",
        ("dashed", "--"): "6,6",
        ("dotted", ":"): "2,2",
        ("dashdot", "-."): "4,4,2,4",
        ("", " ", "None", "none"): None,
    }
)


def get_dasharray(obj):
    """Get an SVG dash array for the given matplotlib linestyle

    Parameters
    ----------
    obj : matplotlib object
        The matplotlib line or path object, which must have a get_linestyle()
        method which returns a valid matplotlib line code

    Returns
    -------
    dasharray : string
        The HTML/SVG dasharray code associated with the object.
    """
    if obj.__dict__.get("_dashSeq", None) is not None:
        return ",".join(map(str, obj._dashSeq))
    else:
        ls = obj.get_linestyle()
        dasharray = LINESTYLES.get(ls, "not found")
        if dasharray == "not found":
            warnings.warn(
                "line style '{0}' not understood: "
                "defaulting to solid line.".format(ls)
            )
            dasharray = LINESTYLES["solid"]
        return dasharray


PATH_DICT = {
    Path.LINETO: "L",
    Path.MOVETO: "M",
    Path.CURVE3: "S",
    Path.CURVE4: "C",
    Path.CLOSEPOLY: "Z",
}


def SVG_path(path, transform=None, simplify=False):
    """Construct the vertices and SVG codes for the path

    Parameters
    ----------
    path : matplotlib.Path object

    transform : matplotlib transform (optional)
        if specified, the path will be transformed before computing the output.

    Returns
    -------
    vertices : array
        The shape (M, 2) array of vertices of the Path. Note that some Path
        codes require multiple vertices, so the length of these vertices may
        be longer than the list of path codes.
    path_codes : list
        A length N list of single-character path codes, N <= M. Each code is
        a single character, in ['L','M','S','C','Z']. See the standard SVG
        path specification for a description of these.
    """
    if transform is not None:
        path = path.transformed(transform)

    vc_tuples = [
        (vertices if path_code != Path.CLOSEPOLY else [], PATH_DICT[path_code])
        for (vertices, path_code) in path.iter_segments(simplify=simplify)
    ]

    if not vc_tuples:
        # empty path is a special case
        return np.zeros((0, 2)), []
    else:
        vertices, codes = zip(*vc_tuples)
        vertices = np.array(list(itertools.chain(*vertices))).reshape(-1, 2)
        return vertices, list(codes)


def get_path_style(path, fill=True):
    """Get the style dictionary for matplotlib path objects"""
    style = {}
    style["alpha"] = path.get_alpha()
    if style["alpha"] is None:
        style["alpha"] = 1
    style["edgecolor"] = export_color(path.get_edgecolor())
    if fill:
        style["facecolor"] = export_color(path.get_facecolor())
    else:
        style["facecolor"] = "none"
    style["edgewidth"] = path.get_linewidth()
    style["dasharray"] = get_dasharray(path)
    style["zorder"] = path.get_zorder()
    return style


def get_line_style(line):
    """Get the style dictionary for matplotlib line objects"""
    style = {}
    style["alpha"] = line.get_alpha()
    if style["alpha"] is None:
        style["alpha"] = 1
    style["color"] = export_color(line.get_color())
    style["linewidth"] = line.get_linewidth()
    style["dasharray"] = get_dasharray(line)
    style["zorder"] = line.get_zorder()
    style["drawstyle"] = line.get_drawstyle()
    return style


def get_marker_style(line):
    """Get the style dictionary for matplotlib marker objects"""
    style = {}
    style["alpha"] = line.get_alpha()
    if style["alpha"] is None:
        style["alpha"] = 1

    style["facecolor"] = export_color(line.get_markerfacecolor())
    style["edgecolor"] = export_color(line.get_markeredgecolor())
    style["edgewidth"] = line.get_markeredgewidth()

    style["marker"] = line.get_marker()
    markerstyle = MarkerStyle(line.get_marker())
    markersize = line.get_markersize()
    markertransform = markerstyle.get_transform() + Affine2D().scale(
        markersize, -markersize
    )
    style["markerpath"] = SVG_path(markerstyle.get_path(), markertransform)
    style["markersize"] = markersize
    style["zorder"] = line.get_zorder()
    return style


def get_text_style(text):
    """Return the text style dict for a text instance"""
    style = {}
    style["alpha"] = text.get_alpha()
    if style["alpha"] is None:
        style["alpha"] = 1
    style["fontsize"] = text.get_size()
    style["color"] = export_color(text.get_color())
    style["halign"] = text.get_horizontalalignment()  # left, center, right
    style["valign"] = text.get_verticalalignment()  # baseline, center, top
    style["malign"] = text._multialignment  # text alignment when '\n' in text
    style["rotation"] = text.get_rotation()
    style["zorder"] = text.get_zorder()
    return style


def get_axis_properties(axis):
    """Return the property dictionary for a matplotlib.Axis instance"""
    props = {}
    label1On = axis._major_tick_kw.get("label1On", True)

    if isinstance(axis, matplotlib.axis.XAxis):
        if label1On:
            props["position"] = "bottom"
        else:
            props["position"] = "top"
    elif isinstance(axis, matplotlib.axis.YAxis):
        if label1On:
            props["position"] = "left"
        else:
            props["position"] = "right"
    else:
        raise ValueError("{0} should be an Axis instance".format(axis))

    # Use tick values if appropriate
    locator = axis.get_major_locator()
    props["nticks"] = len(locator())
    if isinstance(locator, ticker.FixedLocator):
        props["tickvalues"] = list(locator())
    else:
        props["tickvalues"] = None

    # Find tick formats
    formatter = axis.get_major_formatter()
    if isinstance(formatter, ticker.NullFormatter):
        props["tickformat"] = ""
    elif isinstance(formatter, ticker.FixedFormatter):
        props["tickformat"] = list(formatter.seq)
    elif isinstance(formatter, ticker.FuncFormatter):
        props["tickformat"] = list(formatter.func.args[0].values())
    elif not any(label.get_visible() for label in axis.get_ticklabels()):
        props["tickformat"] = ""
    else:
        props["tickformat"] = None

    # Get axis scale
    props["scale"] = axis.get_scale()

    # Get major tick label size (assumes that's all we really care about!)
    labels = axis.get_ticklabels()
    if labels:
        props["fontsize"] = labels[0].get_fontsize()
    else:
        props["fontsize"] = None

    # Get associated grid
    props["grid"] = get_grid_style(axis)

    # get axis visibility
    props["visible"] = axis.get_visible()

    return props


def get_grid_style(axis):
    gridlines = axis.get_gridlines()
    if axis._major_tick_kw["gridOn"] and len(gridlines) > 0:
        color = export_color(gridlines[0].get_color())
        alpha = gridlines[0].get_alpha()
        dasharray = get_dasharray(gridlines[0])
        return dict(gridOn=True, color=color, dasharray=dasharray, alpha=alpha)
    else:
        return {"gridOn": False}


def get_figure_properties(fig):
    return {
        "figwidth": fig.get_figwidth(),
        "figheight": fig.get_figheight(),
        "dpi": fig.dpi,
    }


def get_axes_properties(ax):
    props = {
        "axesbg": export_color(ax.patch.get_facecolor()),
        "axesbgalpha": ax.patch.get_alpha(),
        "bounds": ax.get_position().bounds,
        "dynamic": ax.get_navigate(),
        "axison": ax.axison,
        "frame_on": ax.get_frame_on(),
        "patch_visible": ax.patch.get_visible(),
        "axes": [get_axis_properties(ax.xaxis), get_axis_properties(ax.yaxis)],
    }

    for axname in ["x", "y"]:
        axis = getattr(ax, axname + "axis")
        domain = getattr(ax, "get_{0}lim".format(axname))()
        lim = domain
        if isinstance(axis.converter, matplotlib.dates.DateConverter):
            scale = "date"
            try:
                import pandas as pd
                from pandas.tseries.converter import PeriodConverter
            except ImportError:
                pd = None

            if pd is not None and isinstance(axis.converter, PeriodConverter):
                _dates = [pd.Period(ordinal=int(d), freq=axis.freq) for d in domain]
                domain = [
                    (d.year, d.month - 1, d.day, d.hour, d.minute, d.second, 0)
                    for d in _dates
                ]
            else:
                domain = [
                    (
                        d.year,
                        d.month - 1,
                        d.day,
                        d.hour,
                        d.minute,
                        d.second,
                        d.microsecond * 1e-3,
                    )
                    for d in matplotlib.dates.num2date(domain)
                ]
        else:
            scale = axis.get_scale()

        if scale not in ["date", "linear", "log"]:
            raise ValueError("Unknown axis scale: " "{0}".format(axis.get_scale()))

        props[axname + "scale"] = scale
        props[axname + "lim"] = lim
        props[axname + "domain"] = domain

    return props


def iter_all_children(obj, skipContainers=False):
    """
    Returns an iterator over all childen and nested children using
    obj's get_children() method

    if skipContainers is true, only childless objects are returned.
    """
    if hasattr(obj, "get_children") and len(obj.get_children()) > 0:
        for child in obj.get_children():
            if not skipContainers:
                yield child
            # could use `yield from` in python 3...
            for grandchild in iter_all_children(child, skipContainers):
                yield grandchild
    else:
        yield obj


def get_legend_properties(ax, legend):
    handles, labels = ax.get_legend_handles_labels()
    visible = legend.get_visible()
    return {"handles": handles, "labels": labels, "visible": visible}


def image_to_base64(image):
    """
    Convert a matplotlib image to a base64 png representation

    Parameters
    ----------
    image : matplotlib image object
        The image to be converted.

    Returns
    -------
    image_base64 : string
        The UTF8-encoded base64 string representation of the png image.
    """
    ax = image.axes
    binary_buffer = io.BytesIO()

    # image is saved in axes coordinates: we need to temporarily
    # set the correct limits to get the correct image
    lim = ax.axis()
    ax.axis(image.get_extent())
    image.write_png(binary_buffer)
    ax.axis(lim)

    binary_buffer.seek(0)
    return base64.b64encode(binary_buffer.read()).decode("utf-8")
