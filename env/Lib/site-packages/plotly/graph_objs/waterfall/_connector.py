from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Connector(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "waterfall"
    _path_str = "waterfall.connector"
    _valid_props = {"line", "mode", "visible"}

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.waterfall.connector.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the line color.
                dash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                width
                    Sets the line width (in px).

        Returns
        -------
        plotly.graph_objs.waterfall.connector.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # mode
    # ----
    @property
    def mode(self):
        """
        Sets the shape of connector lines.

        The 'mode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['spanning', 'between']

        Returns
        -------
        Any
        """
        return self["mode"]

    @mode.setter
    def mode(self, val):
        self["mode"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines if connector lines are drawn.

        The 'visible' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.waterfall.connector.Line`
            instance or dict with compatible properties
        mode
            Sets the shape of connector lines.
        visible
            Determines if connector lines are drawn.
        """

    def __init__(self, arg=None, line=None, mode=None, visible=None, **kwargs):
        """
        Construct a new Connector object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.waterfall.Connector`
        line
            :class:`plotly.graph_objects.waterfall.connector.Line`
            instance or dict with compatible properties
        mode
            Sets the shape of connector lines.
        visible
            Determines if connector lines are drawn.

        Returns
        -------
        Connector
        """
        super(Connector, self).__init__("connector")

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
The first argument to the plotly.graph_objs.waterfall.Connector
constructor must be a dict or
an instance of :class:`plotly.graph_objs.waterfall.Connector`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("mode", None)
        _v = mode if mode is not None else _v
        if _v is not None:
            self["mode"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
