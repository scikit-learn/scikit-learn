#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Increasing(_BaseTraceHierarchyType):
    _parent_path_str = "ohlc"
    _path_str = "ohlc.increasing"
    _valid_props = {"line"}

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.ohlc.increasing.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.ohlc.increasing.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.ohlc.increasing.Line`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, line=None, **kwargs):
        """
        Construct a new Increasing object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.ohlc.Increasing`
        line
            :class:`plotly.graph_objects.ohlc.increasing.Line`
            instance or dict with compatible properties

        Returns
        -------
        Increasing
        """
        super().__init__("increasing")
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
The first argument to the plotly.graph_objs.ohlc.Increasing
constructor must be a dict or
an instance of :class:`plotly.graph_objs.ohlc.Increasing`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
