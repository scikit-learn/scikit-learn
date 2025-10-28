#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Delta(_BaseTraceHierarchyType):
    _parent_path_str = "indicator"
    _path_str = "indicator.delta"
    _valid_props = {
        "decreasing",
        "font",
        "increasing",
        "position",
        "prefix",
        "reference",
        "relative",
        "suffix",
        "valueformat",
    }

    @property
    def decreasing(self):
        """
        The 'decreasing' property is an instance of Decreasing
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.delta.Decreasing`
          - A dict of string/value properties that will be passed
            to the Decreasing constructor

        Returns
        -------
        plotly.graph_objs.indicator.delta.Decreasing
        """
        return self["decreasing"]

    @decreasing.setter
    def decreasing(self, val):
        self["decreasing"] = val

    @property
    def font(self):
        """
        Set the font used to display the delta

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.delta.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.indicator.delta.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def increasing(self):
        """
        The 'increasing' property is an instance of Increasing
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.delta.Increasing`
          - A dict of string/value properties that will be passed
            to the Increasing constructor

        Returns
        -------
        plotly.graph_objs.indicator.delta.Increasing
        """
        return self["increasing"]

    @increasing.setter
    def increasing(self, val):
        self["increasing"] = val

    @property
    def position(self):
        """
        Sets the position of delta with respect to the number.

        The 'position' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'bottom', 'left', 'right']

        Returns
        -------
        Any
        """
        return self["position"]

    @position.setter
    def position(self, val):
        self["position"] = val

    @property
    def prefix(self):
        """
        Sets a prefix appearing before the delta.

        The 'prefix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["prefix"]

    @prefix.setter
    def prefix(self, val):
        self["prefix"] = val

    @property
    def reference(self):
        """
        Sets the reference value to compute the delta. By default, it
        is set to the current value.

        The 'reference' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["reference"]

    @reference.setter
    def reference(self, val):
        self["reference"] = val

    @property
    def relative(self):
        """
        Show relative change

        The 'relative' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["relative"]

    @relative.setter
    def relative(self, val):
        self["relative"] = val

    @property
    def suffix(self):
        """
        Sets a suffix appearing next to the delta.

        The 'suffix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["suffix"]

    @suffix.setter
    def suffix(self, val):
        self["suffix"] = val

    @property
    def valueformat(self):
        """
        Sets the value formatting rule using d3 formatting mini-
        languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format.

        The 'valueformat' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["valueformat"]

    @valueformat.setter
    def valueformat(self, val):
        self["valueformat"] = val

    @property
    def _prop_descriptions(self):
        return """\
        decreasing
            :class:`plotly.graph_objects.indicator.delta.Decreasing
            ` instance or dict with compatible properties
        font
            Set the font used to display the delta
        increasing
            :class:`plotly.graph_objects.indicator.delta.Increasing
            ` instance or dict with compatible properties
        position
            Sets the position of delta with respect to the number.
        prefix
            Sets a prefix appearing before the delta.
        reference
            Sets the reference value to compute the delta. By
            default, it is set to the current value.
        relative
            Show relative change
        suffix
            Sets a suffix appearing next to the delta.
        valueformat
            Sets the value formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
        """

    def __init__(
        self,
        arg=None,
        decreasing=None,
        font=None,
        increasing=None,
        position=None,
        prefix=None,
        reference=None,
        relative=None,
        suffix=None,
        valueformat=None,
        **kwargs,
    ):
        """
        Construct a new Delta object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.indicator.Delta`
        decreasing
            :class:`plotly.graph_objects.indicator.delta.Decreasing
            ` instance or dict with compatible properties
        font
            Set the font used to display the delta
        increasing
            :class:`plotly.graph_objects.indicator.delta.Increasing
            ` instance or dict with compatible properties
        position
            Sets the position of delta with respect to the number.
        prefix
            Sets a prefix appearing before the delta.
        reference
            Sets the reference value to compute the delta. By
            default, it is set to the current value.
        relative
            Show relative change
        suffix
            Sets a suffix appearing next to the delta.
        valueformat
            Sets the value formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.

        Returns
        -------
        Delta
        """
        super().__init__("delta")
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
The first argument to the plotly.graph_objs.indicator.Delta
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.Delta`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("decreasing", arg, decreasing)
        self._set_property("font", arg, font)
        self._set_property("increasing", arg, increasing)
        self._set_property("position", arg, position)
        self._set_property("prefix", arg, prefix)
        self._set_property("reference", arg, reference)
        self._set_property("relative", arg, relative)
        self._set_property("suffix", arg, suffix)
        self._set_property("valueformat", arg, valueformat)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
