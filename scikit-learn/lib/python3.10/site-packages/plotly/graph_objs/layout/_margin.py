#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Margin(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.margin"
    _valid_props = {"autoexpand", "b", "l", "pad", "r", "t"}

    @property
    def autoexpand(self):
        """
        Turns on/off margin expansion computations. Legends, colorbars,
        updatemenus, sliders, axis rangeselector and rangeslider are
        allowed to push the margins by defaults.

        The 'autoexpand' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["autoexpand"]

    @autoexpand.setter
    def autoexpand(self, val):
        self["autoexpand"] = val

    @property
    def b(self):
        """
        Sets the bottom margin (in px).

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
        Sets the left margin (in px).

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
    def pad(self):
        """
        Sets the amount of padding (in px) between the plotting area
        and the axis lines

        The 'pad' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    @property
    def r(self):
        """
        Sets the right margin (in px).

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
        Sets the top margin (in px).

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
        autoexpand
            Turns on/off margin expansion computations. Legends,
            colorbars, updatemenus, sliders, axis rangeselector and
            rangeslider are allowed to push the margins by
            defaults.
        b
            Sets the bottom margin (in px).
        l
            Sets the left margin (in px).
        pad
            Sets the amount of padding (in px) between the plotting
            area and the axis lines
        r
            Sets the right margin (in px).
        t
            Sets the top margin (in px).
        """

    def __init__(
        self,
        arg=None,
        autoexpand=None,
        b=None,
        l=None,
        pad=None,
        r=None,
        t=None,
        **kwargs,
    ):
        """
        Construct a new Margin object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Margin`
        autoexpand
            Turns on/off margin expansion computations. Legends,
            colorbars, updatemenus, sliders, axis rangeselector and
            rangeslider are allowed to push the margins by
            defaults.
        b
            Sets the bottom margin (in px).
        l
            Sets the left margin (in px).
        pad
            Sets the amount of padding (in px) between the plotting
            area and the axis lines
        r
            Sets the right margin (in px).
        t
            Sets the top margin (in px).

        Returns
        -------
        Margin
        """
        super().__init__("margin")
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
The first argument to the plotly.graph_objs.layout.Margin
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Margin`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("autoexpand", arg, autoexpand)
        self._set_property("b", arg, b)
        self._set_property("l", arg, l)
        self._set_property("pad", arg, pad)
        self._set_property("r", arg, r)
        self._set_property("t", arg, t)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
