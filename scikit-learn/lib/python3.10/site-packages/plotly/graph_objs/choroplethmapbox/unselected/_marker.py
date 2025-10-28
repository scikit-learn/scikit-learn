#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):
    _parent_path_str = "choroplethmapbox.unselected"
    _path_str = "choroplethmapbox.unselected.marker"
    _valid_props = {"opacity"}

    @property
    def opacity(self):
        """
        Sets the marker opacity of unselected points, applied only when
        a selection exists.

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
    def _prop_descriptions(self):
        return """\
        opacity
            Sets the marker opacity of unselected points, applied
            only when a selection exists.
        """

    def __init__(self, arg=None, opacity=None, **kwargs):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.choroplethmapb
            ox.unselected.Marker`
        opacity
            Sets the marker opacity of unselected points, applied
            only when a selection exists.

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
The first argument to the plotly.graph_objs.choroplethmapbox.unselected.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.choroplethmapbox.unselected.Marker`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("opacity", arg, opacity)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
