#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Bar(_BaseTraceHierarchyType):
    _parent_path_str = "indicator.gauge"
    _path_str = "indicator.gauge.bar"
    _valid_props = {"color", "line", "thickness"}

    @property
    def color(self):
        """
        Sets the background color of the arc.

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
          - An instance of :class:`plotly.graph_objs.indicator.gauge.bar.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.indicator.gauge.bar.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def thickness(self):
        """
        Sets the thickness of the bar as a fraction of the total
        thickness of the gauge.

        The 'thickness' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["thickness"]

    @thickness.setter
    def thickness(self, val):
        self["thickness"] = val

    @property
    def _prop_descriptions(self):
        return """\
        color
            Sets the background color of the arc.
        line
            :class:`plotly.graph_objects.indicator.gauge.bar.Line`
            instance or dict with compatible properties
        thickness
            Sets the thickness of the bar as a fraction of the
            total thickness of the gauge.
        """

    def __init__(self, arg=None, color=None, line=None, thickness=None, **kwargs):
        """
        Construct a new Bar object

        Set the appearance of the gauge's value

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.indicator.gauge.Bar`
        color
            Sets the background color of the arc.
        line
            :class:`plotly.graph_objects.indicator.gauge.bar.Line`
            instance or dict with compatible properties
        thickness
            Sets the thickness of the bar as a fraction of the
            total thickness of the gauge.

        Returns
        -------
        Bar
        """
        super().__init__("bar")
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
The first argument to the plotly.graph_objs.indicator.gauge.Bar
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.gauge.Bar`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("color", arg, color)
        self._set_property("line", arg, line)
        self._set_property("thickness", arg, thickness)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
