#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):
    _parent_path_str = "waterfall.totals"
    _path_str = "waterfall.totals.marker"
    _valid_props = {"color", "line"}

    @property
    def color(self):
        """
        Sets the marker color of all intermediate sums and total
        values.

        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.waterfall.totals.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.waterfall.totals.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def _prop_descriptions(self):
        return """\
        color
            Sets the marker color of all intermediate sums and
            total values.
        line
            :class:`plotly.graph_objects.waterfall.totals.marker.Li
            ne` instance or dict with compatible properties
        """

    def __init__(self, arg=None, color=None, line=None, **kwargs):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.waterfall.totals.Marker`
        color
            Sets the marker color of all intermediate sums and
            total values.
        line
            :class:`plotly.graph_objects.waterfall.totals.marker.Li
            ne` instance or dict with compatible properties

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
The first argument to the plotly.graph_objs.waterfall.totals.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.waterfall.totals.Marker`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("color", arg, color)
        self._set_property("line", arg, line)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
