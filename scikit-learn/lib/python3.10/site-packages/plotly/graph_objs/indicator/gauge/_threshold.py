#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Threshold(_BaseTraceHierarchyType):
    _parent_path_str = "indicator.gauge"
    _path_str = "indicator.gauge.threshold"
    _valid_props = {"line", "thickness", "value"}

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.gauge.threshold.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.indicator.gauge.threshold.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def thickness(self):
        """
        Sets the thickness of the threshold line as a fraction of the
        thickness of the gauge.

        The 'thickness' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["thickness"]

    @thickness.setter
    def thickness(self, val):
        self["thickness"] = val

    @property
    def value(self):
        """
        Sets a treshold value drawn as a line.

        The 'value' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["value"]

    @value.setter
    def value(self, val):
        self["value"] = val

    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.indicator.gauge.threshold.
            Line` instance or dict with compatible properties
        thickness
            Sets the thickness of the threshold line as a fraction
            of the thickness of the gauge.
        value
            Sets a treshold value drawn as a line.
        """

    def __init__(self, arg=None, line=None, thickness=None, value=None, **kwargs):
        """
        Construct a new Threshold object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.indicator.gauge.Threshold`
        line
            :class:`plotly.graph_objects.indicator.gauge.threshold.
            Line` instance or dict with compatible properties
        thickness
            Sets the thickness of the threshold line as a fraction
            of the thickness of the gauge.
        value
            Sets a treshold value drawn as a line.

        Returns
        -------
        Threshold
        """
        super().__init__("threshold")
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
The first argument to the plotly.graph_objs.indicator.gauge.Threshold
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.gauge.Threshold`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._set_property("thickness", arg, thickness)
        self._set_property("value", arg, value)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
