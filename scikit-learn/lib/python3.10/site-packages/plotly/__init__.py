"""
https://plot.ly/python/

Plotly's Python API allows users to programmatically access Plotly's
server resources.

This package is organized as follows:

Subpackages:

- plotly: all functionality that requires access to Plotly's servers

- graph_objs: objects for designing figures and visualizing data

- matplotlylib: tools to convert matplotlib figures

Modules:

- tools: some helpful tools that do not require access to Plotly's servers

- utils: functions that you probably won't need, but that subpackages use

- version: holds the current API version

- exceptions: defines our custom exception classes

"""

from typing import TYPE_CHECKING
from _plotly_utils.importers import relative_import
import importlib.metadata

# This is the version of the plotly package
__version__ = importlib.metadata.version("plotly")
version = __version__

if TYPE_CHECKING:
    from plotly import (
        graph_objs,
        tools,
        utils,
        offline,
        colors,
        io,
        data,
    )
    from plotly.version import __version__

    __all__ = [
        "graph_objs",
        "tools",
        "utils",
        "offline",
        "colors",
        "io",
        "data",
        "__version__",
    ]

    # Set default template (for >= 3.7 this is done in ploty/io/__init__.py)
    from plotly.io import templates

    templates._default = "plotly"
else:
    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [
            ".graph_objs",
            ".graph_objects",
            ".tools",
            ".utils",
            ".offline",
            ".colors",
            ".io",
            ".data",
        ],
        [".version.__version__"],
    )


def plot(data_frame, kind, **kwargs):
    """
    Pandas plotting backend function, not meant to be called directly.
    To activate, set pandas.options.plotting.backend="plotly"
    See https://github.com/pandas-dev/pandas/blob/master/pandas/plotting/__init__.py
    """
    from .express import (
        scatter,
        line,
        area,
        bar,
        box,
        histogram,
        violin,
        strip,
        funnel,
        density_contour,
        density_heatmap,
        imshow,
    )

    if kind == "scatter":
        new_kwargs = {k: kwargs[k] for k in kwargs if k not in ["s", "c"]}
        return scatter(data_frame, **new_kwargs)
    if kind == "line":
        return line(data_frame, **kwargs)
    if kind == "area":
        new_kwargs = {k: kwargs[k] for k in kwargs if k not in ["stacked"]}
        return area(data_frame, **new_kwargs)
    if kind == "bar":
        return bar(data_frame, **kwargs)
    if kind == "barh":
        return bar(data_frame, orientation="h", **kwargs)
    if kind == "box":
        new_kwargs = {k: kwargs[k] for k in kwargs if k not in ["by"]}
        return box(data_frame, **new_kwargs)
    if kind in ["hist", "histogram"]:
        new_kwargs = {k: kwargs[k] for k in kwargs if k not in ["by", "bins"]}
        return histogram(data_frame, **new_kwargs)
    if kind == "violin":
        return violin(data_frame, **kwargs)
    if kind == "strip":
        return strip(data_frame, **kwargs)
    if kind == "funnel":
        return funnel(data_frame, **kwargs)
    if kind == "density_contour":
        return density_contour(data_frame, **kwargs)
    if kind == "density_heatmap":
        return density_heatmap(data_frame, **kwargs)
    if kind == "imshow":
        return imshow(data_frame, **kwargs)
    if kind == "heatmap":
        raise ValueError(
            "kind='heatmap' not supported plotting.backend='plotly'. "
            "Please use kind='imshow' or kind='density_heatmap'."
        )

    raise NotImplementedError(
        "kind='%s' not yet supported for plotting.backend='plotly'" % kind
    )


def boxplot_frame(data_frame, **kwargs):
    """
    Pandas plotting backend function, not meant to be called directly.
    To activate, set pandas.options.plotting.backend="plotly"
    See https://github.com/pandas-dev/pandas/blob/master/pandas/plotting/__init__.py
    """
    from .express import box

    skip = ["by", "column", "ax", "fontsize", "rot", "grid", "figsize", "layout"]
    skip += ["return_type"]
    new_kwargs = {k: kwargs[k] for k in kwargs if k not in skip}
    return box(data_frame, **new_kwargs)


def hist_frame(data_frame, **kwargs):
    """
    Pandas plotting backend function, not meant to be called directly.
    To activate, set pandas.options.plotting.backend="plotly"
    See https://github.com/pandas-dev/pandas/blob/master/pandas/plotting/__init__.py
    """
    from .express import histogram

    skip = ["column", "by", "grid", "xlabelsize", "xrot", "ylabelsize", "yrot"]
    skip += ["ax", "sharex", "sharey", "figsize", "layout", "bins", "legend"]
    new_kwargs = {k: kwargs[k] for k in kwargs if k not in skip}
    return histogram(data_frame, **new_kwargs)


def hist_series(data_frame, **kwargs):
    """
    Pandas plotting backend function, not meant to be called directly.
    To activate, set pandas.options.plotting.backend="plotly"
    See https://github.com/pandas-dev/pandas/blob/master/pandas/plotting/__init__.py
    """
    from .express import histogram

    skip = ["by", "grid", "xlabelsize", "xrot", "ylabelsize", "yrot", "ax"]
    skip += ["figsize", "bins", "legend"]
    new_kwargs = {k: kwargs[k] for k in kwargs if k not in skip}
    return histogram(data_frame, **new_kwargs)


def _jupyter_labextension_paths():
    """Called by Jupyter Lab Server to detect if it is a valid labextension and
    to install the extension.
    """
    return [
        {
            "src": "labextension/static",
            "dest": "jupyterlab-plotly",
        }
    ]
