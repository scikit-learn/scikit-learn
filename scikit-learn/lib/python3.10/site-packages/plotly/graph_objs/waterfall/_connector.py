#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Connector(_BaseTraceHierarchyType):
    _parent_path_str = "waterfall"
    _path_str = "waterfall.connector"
    _valid_props = {"line", "mode", "visible"}

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.waterfall.connector.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.waterfall.connector.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

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
        super().__init__("connector")
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
The first argument to the plotly.graph_objs.waterfall.Connector
constructor must be a dict or
an instance of :class:`plotly.graph_objs.waterfall.Connector`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._set_property("mode", arg, mode)
        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
