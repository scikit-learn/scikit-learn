#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Lonaxis(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.geo"
    _path_str = "layout.geo.lonaxis"
    _valid_props = {
        "dtick",
        "gridcolor",
        "griddash",
        "gridwidth",
        "range",
        "showgrid",
        "tick0",
    }

    @property
    def dtick(self):
        """
        Sets the graticule's longitude/latitude tick step.

        The 'dtick' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["dtick"]

    @dtick.setter
    def dtick(self, val):
        self["dtick"] = val

    @property
    def gridcolor(self):
        """
        Sets the graticule's stroke color.

        The 'gridcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["gridcolor"]

    @gridcolor.setter
    def gridcolor(self, val):
        self["gridcolor"] = val

    @property
    def griddash(self):
        """
        Sets the dash style of lines. Set to a dash type string
        ("solid", "dot", "dash", "longdash", "dashdot", or
        "longdashdot") or a dash length list in px (eg
        "5px,10px,2px,2px").

        The 'griddash' property is an enumeration that may be specified as:
          - One of the following dash styles:
                ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
          - A string containing a dash length list in pixels or percentages
                (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)

        Returns
        -------
        str
        """
        return self["griddash"]

    @griddash.setter
    def griddash(self, val):
        self["griddash"] = val

    @property
    def gridwidth(self):
        """
        Sets the graticule's stroke width (in px).

        The 'gridwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["gridwidth"]

    @gridwidth.setter
    def gridwidth(self, val):
        self["gridwidth"] = val

    @property
    def range(self):
        """
            Sets the range of this axis (in degrees), sets the map's
            clipped coordinates.

            The 'range' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'range[0]' property is a number and may be specified as:
              - An int or float
        (1) The 'range[1]' property is a number and may be specified as:
              - An int or float

            Returns
            -------
            list
        """
        return self["range"]

    @range.setter
    def range(self, val):
        self["range"] = val

    @property
    def showgrid(self):
        """
        Sets whether or not graticule are shown on the map.

        The 'showgrid' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showgrid"]

    @showgrid.setter
    def showgrid(self, val):
        self["showgrid"] = val

    @property
    def tick0(self):
        """
        Sets the graticule's starting tick longitude/latitude.

        The 'tick0' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["tick0"]

    @tick0.setter
    def tick0(self, val):
        self["tick0"] = val

    @property
    def _prop_descriptions(self):
        return """\
        dtick
            Sets the graticule's longitude/latitude tick step.
        gridcolor
            Sets the graticule's stroke color.
        griddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        gridwidth
            Sets the graticule's stroke width (in px).
        range
            Sets the range of this axis (in degrees), sets the
            map's clipped coordinates.
        showgrid
            Sets whether or not graticule are shown on the map.
        tick0
            Sets the graticule's starting tick longitude/latitude.
        """

    def __init__(
        self,
        arg=None,
        dtick=None,
        gridcolor=None,
        griddash=None,
        gridwidth=None,
        range=None,
        showgrid=None,
        tick0=None,
        **kwargs,
    ):
        """
        Construct a new Lonaxis object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.geo.Lonaxis`
        dtick
            Sets the graticule's longitude/latitude tick step.
        gridcolor
            Sets the graticule's stroke color.
        griddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        gridwidth
            Sets the graticule's stroke width (in px).
        range
            Sets the range of this axis (in degrees), sets the
            map's clipped coordinates.
        showgrid
            Sets whether or not graticule are shown on the map.
        tick0
            Sets the graticule's starting tick longitude/latitude.

        Returns
        -------
        Lonaxis
        """
        super().__init__("lonaxis")
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
The first argument to the plotly.graph_objs.layout.geo.Lonaxis
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.geo.Lonaxis`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("dtick", arg, dtick)
        self._set_property("gridcolor", arg, gridcolor)
        self._set_property("griddash", arg, griddash)
        self._set_property("gridwidth", arg, gridwidth)
        self._set_property("range", arg, range)
        self._set_property("showgrid", arg, showgrid)
        self._set_property("tick0", arg, tick0)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
