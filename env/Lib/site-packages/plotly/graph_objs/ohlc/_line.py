from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Line(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "ohlc"
    _path_str = "ohlc.line"
    _valid_props = {"dash", "width"}

    # dash
    # ----
    @property
    def dash(self):
        """
        Sets the dash style of lines. Set to a dash type string
        ("solid", "dot", "dash", "longdash", "dashdot", or
        "longdashdot") or a dash length list in px (eg
        "5px,10px,2px,2px"). Note that this style setting can also be
        set per direction via `increasing.line.dash` and
        `decreasing.line.dash`.

        The 'dash' property is an enumeration that may be specified as:
          - One of the following dash styles:
                ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
          - A string containing a dash length list in pixels or percentages
                (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)

        Returns
        -------
        str
        """
        return self["dash"]

    @dash.setter
    def dash(self, val):
        self["dash"] = val

    # width
    # -----
    @property
    def width(self):
        """
        [object Object] Note that this style setting can also be set
        per direction via `increasing.line.width` and
        `decreasing.line.width`.

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        dash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px"). Note that this style setting can
            also be set per direction via `increasing.line.dash`
            and `decreasing.line.dash`.
        width
            [object Object] Note that this style setting can also
            be set per direction via `increasing.line.width` and
            `decreasing.line.width`.
        """

    def __init__(self, arg=None, dash=None, width=None, **kwargs):
        """
        Construct a new Line object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.ohlc.Line`
        dash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px"). Note that this style setting can
            also be set per direction via `increasing.line.dash`
            and `decreasing.line.dash`.
        width
            [object Object] Note that this style setting can also
            be set per direction via `increasing.line.width` and
            `decreasing.line.width`.

        Returns
        -------
        Line
        """
        super(Line, self).__init__("line")

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
The first argument to the plotly.graph_objs.ohlc.Line
constructor must be a dict or
an instance of :class:`plotly.graph_objs.ohlc.Line`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("dash", None)
        _v = dash if dash is not None else _v
        if _v is not None:
            self["dash"] = _v
        _v = arg.pop("width", None)
        _v = width if width is not None else _v
        if _v is not None:
            self["width"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
