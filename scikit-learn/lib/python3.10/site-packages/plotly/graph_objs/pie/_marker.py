#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):
    _parent_path_str = "pie"
    _path_str = "pie.marker"
    _valid_props = {"colors", "colorssrc", "line", "pattern"}

    @property
    def colors(self):
        """
        Sets the color of each sector. If not specified, the default
        trace color set is used to pick the sector colors.

        The 'colors' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["colors"]

    @colors.setter
    def colors(self, val):
        self["colors"] = val

    @property
    def colorssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `colors`.

        The 'colorssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["colorssrc"]

    @colorssrc.setter
    def colorssrc(self, val):
        self["colorssrc"] = val

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.pie.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.pie.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def pattern(self):
        """
        Sets the pattern within the marker.

        The 'pattern' property is an instance of Pattern
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.pie.marker.Pattern`
          - A dict of string/value properties that will be passed
            to the Pattern constructor

        Returns
        -------
        plotly.graph_objs.pie.marker.Pattern
        """
        return self["pattern"]

    @pattern.setter
    def pattern(self, val):
        self["pattern"] = val

    @property
    def _prop_descriptions(self):
        return """\
        colors
            Sets the color of each sector. If not specified, the
            default trace color set is used to pick the sector
            colors.
        colorssrc
            Sets the source reference on Chart Studio Cloud for
            `colors`.
        line
            :class:`plotly.graph_objects.pie.marker.Line` instance
            or dict with compatible properties
        pattern
            Sets the pattern within the marker.
        """

    def __init__(
        self, arg=None, colors=None, colorssrc=None, line=None, pattern=None, **kwargs
    ):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.pie.Marker`
        colors
            Sets the color of each sector. If not specified, the
            default trace color set is used to pick the sector
            colors.
        colorssrc
            Sets the source reference on Chart Studio Cloud for
            `colors`.
        line
            :class:`plotly.graph_objects.pie.marker.Line` instance
            or dict with compatible properties
        pattern
            Sets the pattern within the marker.

        Returns
        -------
        Marker
        """
        super().__init__("marker")
        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError("""\
The first argument to the plotly.graph_objs.pie.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.pie.Marker`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("colors", arg, colors)
        self._set_property("colorssrc", arg, colorssrc)
        self._set_property("line", arg, line)
        self._set_property("pattern", arg, pattern)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
