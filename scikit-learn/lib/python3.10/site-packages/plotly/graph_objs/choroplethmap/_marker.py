#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):
    _parent_path_str = "choroplethmap"
    _path_str = "choroplethmap.marker"
    _valid_props = {"line", "opacity", "opacitysrc"}

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.choroplethmap.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.choroplethmap.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def opacity(self):
        """
        Sets the opacity of the locations.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    @property
    def opacitysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `opacity`.

        The 'opacitysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["opacitysrc"]

    @opacitysrc.setter
    def opacitysrc(self, val):
        self["opacitysrc"] = val

    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.choroplethmap.marker.Line`
            instance or dict with compatible properties
        opacity
            Sets the opacity of the locations.
        opacitysrc
            Sets the source reference on Chart Studio Cloud for
            `opacity`.
        """

    def __init__(self, arg=None, line=None, opacity=None, opacitysrc=None, **kwargs):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.choroplethmap.Marker`
        line
            :class:`plotly.graph_objects.choroplethmap.marker.Line`
            instance or dict with compatible properties
        opacity
            Sets the opacity of the locations.
        opacitysrc
            Sets the source reference on Chart Studio Cloud for
            `opacity`.

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
The first argument to the plotly.graph_objs.choroplethmap.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.choroplethmap.Marker`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._set_property("opacity", arg, opacity)
        self._set_property("opacitysrc", arg, opacitysrc)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
