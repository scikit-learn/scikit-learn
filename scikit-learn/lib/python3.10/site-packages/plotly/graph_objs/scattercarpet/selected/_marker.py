#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):
    _parent_path_str = "scattercarpet.selected"
    _path_str = "scattercarpet.selected.marker"
    _valid_props = {"color", "opacity", "size"}

    @property
    def color(self):
        """
        Sets the marker color of selected points.

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
    def opacity(self):
        """
        Sets the marker opacity of selected points.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    @property
    def size(self):
        """
        Sets the marker size of selected points.

        The 'size' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    @property
    def _prop_descriptions(self):
        return """\
        color
            Sets the marker color of selected points.
        opacity
            Sets the marker opacity of selected points.
        size
            Sets the marker size of selected points.
        """

    def __init__(self, arg=None, color=None, opacity=None, size=None, **kwargs):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.scattercarpet.
            selected.Marker`
        color
            Sets the marker color of selected points.
        opacity
            Sets the marker opacity of selected points.
        size
            Sets the marker size of selected points.

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
The first argument to the plotly.graph_objs.scattercarpet.selected.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scattercarpet.selected.Marker`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("color", arg, color)
        self._set_property("opacity", arg, opacity)
        self._set_property("size", arg, size)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
