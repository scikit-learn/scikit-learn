from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Margin(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.margin"
    _valid_props = {"autoexpand", "b", "l", "pad", "r", "t"}

    # autoexpand
    # ----------
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

    # b
    # -
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

    # l
    # -
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

    # pad
    # ---
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

    # r
    # -
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

    # t
    # -
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

    # Self properties description
    # ---------------------------
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
        super(Margin, self).__init__("margin")

        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        # Validate arg
        # ------------
        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError(
                """\
The first argument to the plotly.graph_objs.layout.Margin
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Margin`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("autoexpand", None)
        _v = autoexpand if autoexpand is not None else _v
        if _v is not None:
            self["autoexpand"] = _v
        _v = arg.pop("b", None)
        _v = b if b is not None else _v
        if _v is not None:
            self["b"] = _v
        _v = arg.pop("l", None)
        _v = l if l is not None else _v
        if _v is not None:
            self["l"] = _v
        _v = arg.pop("pad", None)
        _v = pad if pad is not None else _v
        if _v is not None:
            self["pad"] = _v
        _v = arg.pop("r", None)
        _v = r if r is not None else _v
        if _v is not None:
            self["r"] = _v
        _v = arg.pop("t", None)
        _v = t if t is not None else _v
        if _v is not None:
            self["t"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
