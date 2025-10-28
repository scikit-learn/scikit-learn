#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Line(_BaseTraceHierarchyType):
    _parent_path_str = "candlestick"
    _path_str = "candlestick.line"
    _valid_props = {"width"}

    @property
    def width(self):
        """
        Sets the width (in px) of line bounding the box(es). Note that
        this style setting can also be set per direction via
        `increasing.line.width` and `decreasing.line.width`.

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

    @property
    def _prop_descriptions(self):
        return """\
        width
            Sets the width (in px) of line bounding the box(es).
            Note that this style setting can also be set per
            direction via `increasing.line.width` and
            `decreasing.line.width`.
        """

    def __init__(self, arg=None, width=None, **kwargs):
        """
        Construct a new Line object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.candlestick.Line`
        width
            Sets the width (in px) of line bounding the box(es).
            Note that this style setting can also be set per
            direction via `increasing.line.width` and
            `decreasing.line.width`.

        Returns
        -------
        Line
        """
        super().__init__("line")
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
The first argument to the plotly.graph_objs.candlestick.Line
constructor must be a dict or
an instance of :class:`plotly.graph_objs.candlestick.Line`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("width", arg, width)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
