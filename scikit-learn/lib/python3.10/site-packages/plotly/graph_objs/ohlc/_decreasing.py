#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Decreasing(_BaseTraceHierarchyType):
    _parent_path_str = "ohlc"
    _path_str = "ohlc.decreasing"
    _valid_props = {"line"}

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.ohlc.decreasing.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.ohlc.decreasing.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.ohlc.decreasing.Line`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, line=None, **kwargs):
        """
        Construct a new Decreasing object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.ohlc.Decreasing`
        line
            :class:`plotly.graph_objects.ohlc.decreasing.Line`
            instance or dict with compatible properties

        Returns
        -------
        Decreasing
        """
        super().__init__("decreasing")
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
The first argument to the plotly.graph_objs.ohlc.Decreasing
constructor must be a dict or
an instance of :class:`plotly.graph_objs.ohlc.Decreasing`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
