#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Pad(_BaseTraceHierarchyType):
    _parent_path_str = "treemap.marker"
    _path_str = "treemap.marker.pad"
    _valid_props = {"b", "l", "r", "t"}

    @property
    def b(self):
        """
        Sets the padding form the bottom (in px).

        The 'b' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["b"]

    @b.setter
    def b(self, val):
        self["b"] = val

    @property
    def l(self):
        """
        Sets the padding form the left (in px).

        The 'l' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["l"]

    @l.setter
    def l(self, val):
        self["l"] = val

    @property
    def r(self):
        """
        Sets the padding form the right (in px).

        The 'r' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["r"]

    @r.setter
    def r(self, val):
        self["r"] = val

    @property
    def t(self):
        """
        Sets the padding form the top (in px).

        The 't' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["t"]

    @t.setter
    def t(self, val):
        self["t"] = val

    @property
    def _prop_descriptions(self):
        return """\
        b
            Sets the padding form the bottom (in px).
        l
            Sets the padding form the left (in px).
        r
            Sets the padding form the right (in px).
        t
            Sets the padding form the top (in px).
        """

    def __init__(self, arg=None, b=None, l=None, r=None, t=None, **kwargs):
        """
        Construct a new Pad object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.treemap.marker.Pad`
        b
            Sets the padding form the bottom (in px).
        l
            Sets the padding form the left (in px).
        r
            Sets the padding form the right (in px).
        t
            Sets the padding form the top (in px).

        Returns
        -------
        Pad
        """
        super().__init__("pad")
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
The first argument to the plotly.graph_objs.treemap.marker.Pad
constructor must be a dict or
an instance of :class:`plotly.graph_objs.treemap.marker.Pad`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("b", arg, b)
        self._set_property("l", arg, l)
        self._set_property("r", arg, r)
        self._set_property("t", arg, t)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
