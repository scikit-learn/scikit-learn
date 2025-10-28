import warnings
import itertools
from contextlib import contextmanager
from packaging.version import Version

import numpy as np
import matplotlib as mpl
from matplotlib import transforms

from .. import utils


class Renderer(object):
    @staticmethod
    def ax_zoomable(ax):
        return bool(ax and ax.get_navigate())

    @staticmethod
    def ax_has_xgrid(ax):
        return bool(ax and ax.xaxis._gridOnMajor and ax.yaxis.get_gridlines())

    @staticmethod
    def ax_has_ygrid(ax):
        return bool(ax and ax.yaxis._gridOnMajor and ax.yaxis.get_gridlines())

    @property
    def current_ax_zoomable(self):
        return self.ax_zoomable(self._current_ax)

    @property
    def current_ax_has_xgrid(self):
        return self.ax_has_xgrid(self._current_ax)

    @property
    def current_ax_has_ygrid(self):
        return self.ax_has_ygrid(self._current_ax)

    @contextmanager
    def draw_figure(self, fig, props):
        if hasattr(self, "_current_fig") and self._current_fig is not None:
            warnings.warn("figure embedded in figure: something is wrong")
        self._current_fig = fig
        self._fig_props = props
        self.open_figure(fig=fig, props=props)
        yield
        self.close_figure(fig=fig)
        self._current_fig = None
        self._fig_props = {}

    @contextmanager
    def draw_axes(self, ax, props):
        if hasattr(self, "_current_ax") and self._current_ax is not None:
            warnings.warn("axes embedded in axes: something is wrong")
        self._current_ax = ax
        self._ax_props = props
        self.open_axes(ax=ax, props=props)
        yield
        self.close_axes(ax=ax)
        self._current_ax = None
        self._ax_props = {}

    @contextmanager
    def draw_legend(self, legend, props):
        self._current_legend = legend
        self._legend_props = props
        self.open_legend(legend=legend, props=props)
        yield
        self.close_legend(legend=legend)
        self._current_legend = None
        self._legend_props = {}

    # Following are the functions which should be overloaded in subclasses

    def open_figure(self, fig, props):
        """
        Begin commands for a particular figure.

        Parameters
        ----------
        fig : matplotlib.Figure
            The Figure which will contain the ensuing axes and elements
        props : dictionary
            The dictionary of figure properties
        """
        pass

    def close_figure(self, fig):
        """
        Finish commands for a particular figure.

        Parameters
        ----------
        fig : matplotlib.Figure
            The figure which is finished being drawn.
        """
        pass

    def open_axes(self, ax, props):
        """
        Begin commands for a particular axes.

        Parameters
        ----------
        ax : matplotlib.Axes
            The Axes which will contain the ensuing axes and elements
        props : dictionary
            The dictionary of axes properties
        """
        pass

    def close_axes(self, ax):
        """
        Finish commands for a particular axes.

        Parameters
        ----------
        ax : matplotlib.Axes
            The Axes which is finished being drawn.
        """
        pass

    def open_legend(self, legend, props):
        """
        Beging commands for a particular legend.

        Parameters
        ----------
        legend : matplotlib.legend.Legend
                The Legend that will contain the ensuing elements
        props : dictionary
                The dictionary of legend properties
        """
        pass

    def close_legend(self, legend):
        """
        Finish commands for a particular legend.

        Parameters
        ----------
        legend : matplotlib.legend.Legend
                The Legend which is finished being drawn
        """
        pass

    def draw_marked_line(
        self, data, coordinates, linestyle, markerstyle, label, mplobj=None
    ):
        """Draw a line that also has markers.

        If this isn't reimplemented by a renderer object, by default, it will
        make a call to BOTH draw_line and draw_markers when both markerstyle
        and linestyle are not None in the same Line2D object.

        """
        if linestyle is not None:
            self.draw_line(data, coordinates, linestyle, label, mplobj)
        if markerstyle is not None:
            self.draw_markers(data, coordinates, markerstyle, label, mplobj)

    def draw_line(self, data, coordinates, style, label, mplobj=None):
        """
        Draw a line. By default, draw the line via the draw_path() command.
        Some renderers might wish to override this and provide more
        fine-grained behavior.

        In matplotlib, lines are generally created via the plt.plot() command,
        though this command also can create marker collections.

        Parameters
        ----------
        data : array_like
            A shape (N, 2) array of datapoints.
        coordinates : string
            A string code, which should be either 'data' for data coordinates,
            or 'figure' for figure (pixel) coordinates.
        style : dictionary
            a dictionary specifying the appearance of the line.
        mplobj : matplotlib object
            the matplotlib plot element which generated this line
        """
        pathcodes = ["M"] + (data.shape[0] - 1) * ["L"]
        pathstyle = dict(facecolor="none", **style)
        pathstyle["edgecolor"] = pathstyle.pop("color")
        pathstyle["edgewidth"] = pathstyle.pop("linewidth")
        self.draw_path(
            data=data,
            coordinates=coordinates,
            pathcodes=pathcodes,
            style=pathstyle,
            mplobj=mplobj,
        )

    @staticmethod
    def _iter_path_collection(paths, path_transforms, offsets, styles):
        """Build an iterator over the elements of the path collection"""
        N = max(len(paths), len(offsets))

        # Before mpl 1.4.0, path_transform can be a false-y value, not a valid
        # transformation matrix.
        if Version(mpl.__version__) < Version("1.4.0"):
            if path_transforms is None:
                path_transforms = [np.eye(3)]

        edgecolor = styles["edgecolor"]
        if np.size(edgecolor) == 0:
            edgecolor = ["none"]
        facecolor = styles["facecolor"]
        if np.size(facecolor) == 0:
            facecolor = ["none"]

        elements = [
            paths,
            path_transforms,
            offsets,
            edgecolor,
            styles["linewidth"],
            facecolor,
        ]

        it = itertools
        return it.islice(zip(*map(it.cycle, elements)), N)

    def draw_path_collection(
        self,
        paths,
        path_coordinates,
        path_transforms,
        offsets,
        offset_coordinates,
        offset_order,
        styles,
        mplobj=None,
    ):
        """
        Draw a collection of paths. The paths, offsets, and styles are all
        iterables, and the number of paths is max(len(paths), len(offsets)).

        By default, this is implemented via multiple calls to the draw_path()
        function. For efficiency, Renderers may choose to customize this
        implementation.

        Examples of path collections created by matplotlib are scatter plots,
        histograms, contour plots, and many others.

        Parameters
        ----------
        paths : list
            list of tuples, where each tuple has two elements:
            (data, pathcodes).  See draw_path() for a description of these.
        path_coordinates: string
            the coordinates code for the paths, which should be either
            'data' for data coordinates, or 'figure' for figure (pixel)
            coordinates.
        path_transforms: array_like
            an array of shape (*, 3, 3), giving a series of 2D Affine
            transforms for the paths. These encode translations, rotations,
            and scalings in the standard way.
        offsets: array_like
            An array of offsets of shape (N, 2)
        offset_coordinates : string
            the coordinates code for the offsets, which should be either
            'data' for data coordinates, or 'figure' for figure (pixel)
            coordinates.
        offset_order : string
            either "before" or "after". This specifies whether the offset
            is applied before the path transform, or after.  The matplotlib
            backend equivalent is "before"->"data", "after"->"screen".
        styles: dictionary
            A dictionary in which each value is a list of length N, containing
            the style(s) for the paths.
        mplobj : matplotlib object
            the matplotlib plot element which generated this collection
        """
        if offset_order == "before":
            raise NotImplementedError("offset before transform")

        for tup in self._iter_path_collection(paths, path_transforms, offsets, styles):
            (path, path_transform, offset, ec, lw, fc) = tup
            vertices, pathcodes = path
            path_transform = transforms.Affine2D(path_transform)
            vertices = path_transform.transform(vertices)
            # This is a hack:
            if path_coordinates == "figure":
                path_coordinates = "points"
            style = {
                "edgecolor": utils.export_color(ec),
                "facecolor": utils.export_color(fc),
                "edgewidth": lw,
                "dasharray": "10,0",
                "alpha": styles["alpha"],
                "zorder": styles["zorder"],
            }
            self.draw_path(
                data=vertices,
                coordinates=path_coordinates,
                pathcodes=pathcodes,
                style=style,
                offset=offset,
                offset_coordinates=offset_coordinates,
                mplobj=mplobj,
            )

    def draw_markers(self, data, coordinates, style, label, mplobj=None):
        """
        Draw a set of markers. By default, this is done by repeatedly
        calling draw_path(), but renderers should generally overload
        this method to provide a more efficient implementation.

        In matplotlib, markers are created using the plt.plot() command.

        Parameters
        ----------
        data : array_like
            A shape (N, 2) array of datapoints.
        coordinates : string
            A string code, which should be either 'data' for data coordinates,
            or 'figure' for figure (pixel) coordinates.
        style : dictionary
            a dictionary specifying the appearance of the markers.
        mplobj : matplotlib object
            the matplotlib plot element which generated this marker collection
        """
        vertices, pathcodes = style["markerpath"]
        pathstyle = dict(
            (key, style[key])
            for key in ["alpha", "edgecolor", "facecolor", "zorder", "edgewidth"]
        )
        pathstyle["dasharray"] = "10,0"
        for vertex in data:
            self.draw_path(
                data=vertices,
                coordinates="points",
                pathcodes=pathcodes,
                style=pathstyle,
                offset=vertex,
                offset_coordinates=coordinates,
                mplobj=mplobj,
            )

    def draw_text(
        self, text, position, coordinates, style, text_type=None, mplobj=None
    ):
        """
        Draw text on the image.

        Parameters
        ----------
        text : string
            The text to draw
        position : tuple
            The (x, y) position of the text
        coordinates : string
            A string code, which should be either 'data' for data coordinates,
            or 'figure' for figure (pixel) coordinates.
        style : dictionary
            a dictionary specifying the appearance of the text.
        text_type : string or None
            if specified, a type of text such as "xlabel", "ylabel", "title"
        mplobj : matplotlib object
            the matplotlib plot element which generated this text
        """
        raise NotImplementedError()

    def draw_path(
        self,
        data,
        coordinates,
        pathcodes,
        style,
        offset=None,
        offset_coordinates="data",
        mplobj=None,
    ):
        """
        Draw a path.

        In matplotlib, paths are created by filled regions, histograms,
        contour plots, patches, etc.

        Parameters
        ----------
        data : array_like
            A shape (N, 2) array of datapoints.
        coordinates : string
            A string code, which should be either 'data' for data coordinates,
            'figure' for figure (pixel) coordinates, or "points" for raw
            point coordinates (useful in conjunction with offsets, below).
        pathcodes : list
            A list of single-character SVG pathcodes associated with the data.
            Path codes are one of ['M', 'm', 'L', 'l', 'Q', 'q', 'T', 't',
                                   'S', 's', 'C', 'c', 'Z', 'z']
            See the SVG specification for details.  Note that some path codes
            consume more than one datapoint (while 'Z' consumes none), so
            in general, the length of the pathcodes list will not be the same
            as that of the data array.
        style : dictionary
            a dictionary specifying the appearance of the line.
        offset : list (optional)
            the (x, y) offset of the path. If not given, no offset will
            be used.
        offset_coordinates : string (optional)
            A string code, which should be either 'data' for data coordinates,
            or 'figure' for figure (pixel) coordinates.
        mplobj : matplotlib object
            the matplotlib plot element which generated this path
        """
        raise NotImplementedError()

    def draw_image(self, imdata, extent, coordinates, style, mplobj=None):
        """
        Draw an image.

        Parameters
        ----------
        imdata : string
            base64 encoded png representation of the image
        extent : list
            the axes extent of the image: [xmin, xmax, ymin, ymax]
        coordinates: string
            A string code, which should be either 'data' for data coordinates,
            or 'figure' for figure (pixel) coordinates.
        style : dictionary
            a dictionary specifying the appearance of the image
        mplobj : matplotlib object
            the matplotlib plot object which generated this image
        """
        raise NotImplementedError()
