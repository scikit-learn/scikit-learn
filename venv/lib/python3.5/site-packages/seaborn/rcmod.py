"""Functions that alter the matplotlib rc dictionary on the fly."""
from distutils.version import LooseVersion
import functools

import numpy as np
import matplotlib as mpl

from . import palettes, _orig_rc_params


mpl_ge_150 = LooseVersion(mpl.__version__) >= '1.5.0'


__all__ = ["set", "reset_defaults", "reset_orig",
           "axes_style", "set_style", "plotting_context", "set_context",
           "set_palette"]


_style_keys = (

    "axes.facecolor",
    "axes.edgecolor",
    "axes.grid",
    "axes.axisbelow",
    "axes.linewidth",
    "axes.labelcolor",

    "figure.facecolor",

    "grid.color",
    "grid.linestyle",

    "text.color",

    "xtick.color",
    "ytick.color",
    "xtick.direction",
    "ytick.direction",
    "xtick.major.size",
    "ytick.major.size",
    "xtick.minor.size",
    "ytick.minor.size",

    "legend.frameon",
    "legend.numpoints",
    "legend.scatterpoints",

    "lines.solid_capstyle",

    "image.cmap",
    "font.family",
    "font.sans-serif",
    )

_context_keys = (

    "font.size",
    "axes.labelsize",
    "axes.titlesize",
    "xtick.labelsize",
    "ytick.labelsize",
    "legend.fontsize",

    "grid.linewidth",
    "lines.linewidth",
    "patch.linewidth",
    "lines.markersize",
    "lines.markeredgewidth",

    "xtick.major.width",
    "ytick.major.width",
    "xtick.minor.width",
    "ytick.minor.width",

    "xtick.major.pad",
    "ytick.major.pad"
    )


def set(context="notebook", style="darkgrid", palette="deep",
        font="sans-serif", font_scale=1, color_codes=False, rc=None):
    """Set aesthetic parameters in one step.

    Each set of parameters can be set directly or temporarily, see the
    referenced functions below for more information.

    Parameters
    ----------
    context : string or dict
        Plotting context parameters, see :func:`plotting_context`
    style : string or dict
        Axes style parameters, see :func:`axes_style`
    palette : string or sequence
        Color palette, see :func:`color_palette`
    font : string
        Font family, see matplotlib font manager.
    font_scale : float, optional
        Separate scaling factor to independently scale the size of the
        font elements.
    color_codes : bool
        If ``True`` and ``palette`` is a seaborn palette, remap the shorthand
        color codes (e.g. "b", "g", "r", etc.) to the colors from this palette.
    rc : dict or None
        Dictionary of rc parameter mappings to override the above.

    """
    set_context(context, font_scale)
    set_style(style, rc={"font.family": font})
    set_palette(palette, color_codes=color_codes)
    if rc is not None:
        mpl.rcParams.update(rc)


def reset_defaults():
    """Restore all RC params to default settings."""
    mpl.rcParams.update(mpl.rcParamsDefault)


def reset_orig():
    """Restore all RC params to original settings (respects custom rc)."""
    mpl.rcParams.update(_orig_rc_params)


def axes_style(style=None, rc=None):
    """Return a parameter dict for the aesthetic style of the plots.

    This affects things like the color of the axes, whether a grid is
    enabled by default, and other aesthetic elements.

    This function returns an object that can be used in a ``with`` statement
    to temporarily change the style parameters.

    Parameters
    ----------
    style : dict, None, or one of {darkgrid, whitegrid, dark, white, ticks}
        A dictionary of parameters or the name of a preconfigured set.
    rc : dict, optional
        Parameter mappings to override the values in the preset seaborn
        style dictionaries. This only updates parameters that are
        considered part of the style definition.

    Examples
    --------
    >>> st = axes_style("whitegrid")

    >>> set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})

    >>> import matplotlib.pyplot as plt
    >>> with axes_style("white"):
    ...     f, ax = plt.subplots()
    ...     ax.plot(x, y)               # doctest: +SKIP

    See Also
    --------
    set_style : set the matplotlib parameters for a seaborn theme
    plotting_context : return a parameter dict to to scale plot elements
    color_palette : define the color palette for a plot

    """
    if style is None:
        style_dict = {k: mpl.rcParams[k] for k in _style_keys}

    elif isinstance(style, dict):
        style_dict = style

    else:
        styles = ["white", "dark", "whitegrid", "darkgrid", "ticks"]
        if style not in styles:
            raise ValueError("style must be one of %s" % ", ".join(styles))

        # Define colors here
        dark_gray = ".15"
        light_gray = ".8"

        # Common parameters
        style_dict = {
            "figure.facecolor": "white",
            "text.color": dark_gray,
            "axes.labelcolor": dark_gray,
            "legend.frameon": False,
            "legend.numpoints": 1,
            "legend.scatterpoints": 1,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.color": dark_gray,
            "ytick.color": dark_gray,
            "axes.axisbelow": True,
            "image.cmap": "rocket",
            "font.family": ["sans-serif"],
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans",
                                "Bitstream Vera Sans", "sans-serif"],
            "grid.linestyle": "-",
            "lines.solid_capstyle": "round",
            }

        # Set grid on or off
        if "grid" in style:
            style_dict.update({
                "axes.grid": True,
                })
        else:
            style_dict.update({
                "axes.grid": False,
                })

        # Set the color of the background, spines, and grids
        if style.startswith("dark"):
            style_dict.update({
                "axes.facecolor": "#EAEAF2",
                "axes.edgecolor": "white",
                "axes.linewidth": 0,
                "grid.color": "white",
                })

        elif style == "whitegrid":
            style_dict.update({
                "axes.facecolor": "white",
                "axes.edgecolor": light_gray,
                "axes.linewidth": 1,
                "grid.color": light_gray,
                })

        elif style in ["white", "ticks"]:
            style_dict.update({
                "axes.facecolor": "white",
                "axes.edgecolor": dark_gray,
                "axes.linewidth": 1.25,
                "grid.color": light_gray,
                })

        # Show or hide the axes ticks
        if style == "ticks":
            style_dict.update({
                "xtick.major.size": 6,
                "ytick.major.size": 6,
                "xtick.minor.size": 3,
                "ytick.minor.size": 3,
                })
        else:
            style_dict.update({
                "xtick.major.size": 0,
                "ytick.major.size": 0,
                "xtick.minor.size": 0,
                "ytick.minor.size": 0,
                })

    # Override these settings with the provided rc dictionary
    if rc is not None:
        rc = {k: v for k, v in rc.items() if k in _style_keys}
        style_dict.update(rc)

    # Wrap in an _AxesStyle object so this can be used in a with statement
    style_object = _AxesStyle(style_dict)

    return style_object


def set_style(style=None, rc=None):
    """Set the aesthetic style of the plots.

    This affects things like the color of the axes, whether a grid is
    enabled by default, and other aesthetic elements.

    Parameters
    ----------
    style : dict, None, or one of {darkgrid, whitegrid, dark, white, ticks}
        A dictionary of parameters or the name of a preconfigured set.
    rc : dict, optional
        Parameter mappings to override the values in the preset seaborn
        style dictionaries. This only updates parameters that are
        considered part of the style definition.

    Examples
    --------
    >>> set_style("whitegrid")

    >>> set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})

    See Also
    --------
    axes_style : return a dict of parameters or use in a ``with`` statement
                 to temporarily set the style.
    set_context : set parameters to scale plot elements
    set_palette : set the default color palette for figures

    """
    style_object = axes_style(style, rc)
    mpl.rcParams.update(style_object)


def plotting_context(context=None, font_scale=1, rc=None):
    """Return a parameter dict to scale elements of the figure.

    This affects things like the size of the labels, lines, and other
    elements of the plot, but not the overall style. The base context
    is "notebook", and the other contexts are "paper", "talk", and "poster",
    which are version of the notebook parameters scaled by .8, 1.3, and 1.6,
    respectively.

    This function returns an object that can be used in a ``with`` statement
    to temporarily change the context parameters.

    Parameters
    ----------
    context : dict, None, or one of {paper, notebook, talk, poster}
        A dictionary of parameters or the name of a preconfigured set.
    font_scale : float, optional
        Separate scaling factor to independently scale the size of the
        font elements.
    rc : dict, optional
        Parameter mappings to override the values in the preset seaborn
        context dictionaries. This only updates parameters that are
        considered part of the context definition.

    Examples
    --------
    >>> c = plotting_context("poster")

    >>> c = plotting_context("notebook", font_scale=1.5)

    >>> c = plotting_context("talk", rc={"lines.linewidth": 2})

    >>> import matplotlib.pyplot as plt
    >>> with plotting_context("paper"):
    ...     f, ax = plt.subplots()
    ...     ax.plot(x, y)                 # doctest: +SKIP

    See Also
    --------
    set_context : set the matplotlib parameters to scale plot elements
    axes_style : return a dict of parameters defining a figure style
    color_palette : define the color palette for a plot

    """
    if context is None:
        context_dict = {k: mpl.rcParams[k] for k in _context_keys}

    elif isinstance(context, dict):
        context_dict = context

    else:

        contexts = ["paper", "notebook", "talk", "poster"]
        if context not in contexts:
            raise ValueError("context must be in %s" % ", ".join(contexts))

        # Set up dictionary of default parameters
        base_context = {

            "font.size": 12,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,

            "grid.linewidth": 1,
            "lines.linewidth": 1.75,
            "patch.linewidth": .3,
            "lines.markersize": 7,
            "lines.markeredgewidth": 0,

            "xtick.major.width": 1,
            "ytick.major.width": 1,
            "xtick.minor.width": .5,
            "ytick.minor.width": .5,

            "xtick.major.pad": 7,
            "ytick.major.pad": 7,

            }

        # Scale all the parameters by the same factor depending on the context
        scaling = dict(paper=.8, notebook=1, talk=1.3, poster=1.6)[context]
        context_dict = {k: v * scaling for k, v in base_context.items()}

        # Now independently scale the fonts
        font_keys = ["axes.labelsize", "axes.titlesize", "legend.fontsize",
                     "xtick.labelsize", "ytick.labelsize", "font.size"]
        font_dict = {k: context_dict[k] * font_scale for k in font_keys}
        context_dict.update(font_dict)

    # Implement hack workaround for matplotlib bug
    # See https://github.com/mwaskom/seaborn/issues/344
    # There is a bug in matplotlib 1.4.2 that makes points invisible when
    # they don't have an edgewidth. It will supposedly be fixed in 1.4.3.
    if mpl.__version__ == "1.4.2":
        context_dict["lines.markeredgewidth"] = 0.01

    # Override these settings with the provided rc dictionary
    if rc is not None:
        rc = {k: v for k, v in rc.items() if k in _context_keys}
        context_dict.update(rc)

    # Wrap in a _PlottingContext object so this can be used in a with statement
    context_object = _PlottingContext(context_dict)

    return context_object


def set_context(context=None, font_scale=1, rc=None):
    """Set the plotting context parameters.

    This affects things like the size of the labels, lines, and other
    elements of the plot, but not the overall style. The base context
    is "notebook", and the other contexts are "paper", "talk", and "poster",
    which are version of the notebook parameters scaled by .8, 1.3, and 1.6,
    respectively.

    Parameters
    ----------
    context : dict, None, or one of {paper, notebook, talk, poster}
        A dictionary of parameters or the name of a preconfigured set.
    font_scale : float, optional
        Separate scaling factor to independently scale the size of the
        font elements.
    rc : dict, optional
        Parameter mappings to override the values in the preset seaborn
        context dictionaries. This only updates parameters that are
        considered part of the context definition.

    Examples
    --------
    >>> set_context("paper")

    >>> set_context("talk", font_scale=1.4)

    >>> set_context("talk", rc={"lines.linewidth": 2})

    See Also
    --------
    plotting_context : return a dictionary of rc parameters, or use in
                       a ``with`` statement to temporarily set the context.
    set_style : set the default parameters for figure style
    set_palette : set the default color palette for figures

    """
    context_object = plotting_context(context, font_scale, rc)
    mpl.rcParams.update(context_object)


class _RCAesthetics(dict):
    def __enter__(self):
        rc = mpl.rcParams
        self._orig = {k: rc[k] for k in self._keys}
        self._set(self)

    def __exit__(self, exc_type, exc_value, exc_tb):
        self._set(self._orig)

    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with self:
                return func(*args, **kwargs)
        return wrapper


class _AxesStyle(_RCAesthetics):
    """Light wrapper on a dict to set style temporarily."""
    _keys = _style_keys
    _set = staticmethod(set_style)


class _PlottingContext(_RCAesthetics):
    """Light wrapper on a dict to set context temporarily."""
    _keys = _context_keys
    _set = staticmethod(set_context)


def set_palette(palette, n_colors=None, desat=None, color_codes=False):
    """Set the matplotlib color cycle using a seaborn palette.

    Parameters
    ----------
    palette : hls | husl | matplotlib colormap | seaborn color palette
        Palette definition. Should be something that :func:`color_palette`
        can process.
    n_colors : int
        Number of colors in the cycle. The default number of colors will depend
        on the format of ``palette``, see the :func:`color_palette`
        documentation for more information.
    desat : float
        Proportion to desaturate each color by.
    color_codes : bool
        If ``True`` and ``palette`` is a seaborn palette, remap the shorthand
        color codes (e.g. "b", "g", "r", etc.) to the colors from this palette.

    Examples
    --------
    >>> set_palette("Reds")

    >>> set_palette("Set1", 8, .75)

    See Also
    --------
    color_palette : build a color palette or set the color cycle temporarily
                    in a ``with`` statement.
    set_context : set parameters to scale plot elements
    set_style : set the default parameters for figure style

    """
    colors = palettes.color_palette(palette, n_colors, desat)
    if mpl_ge_150:
        from cycler import cycler
        cyl = cycler('color', colors)
        mpl.rcParams['axes.prop_cycle'] = cyl
    else:
        mpl.rcParams["axes.color_cycle"] = list(colors)
    mpl.rcParams["patch.facecolor"] = colors[0]
    if color_codes:
        palettes.set_color_codes(palette)
