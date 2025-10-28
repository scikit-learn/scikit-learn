#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Number(_BaseTraceHierarchyType):
    _parent_path_str = "indicator"
    _path_str = "indicator.number"
    _valid_props = {"font", "prefix", "suffix", "valueformat"}

    @property
    def font(self):
        """
        Set the font used to display main number

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.number.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.indicator.number.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def prefix(self):
        """
        Sets a prefix appearing before the number.

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
    def suffix(self):
        """
        Sets a suffix appearing next to the number.

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
        font
            Set the font used to display main number
        prefix
            Sets a prefix appearing before the number.
        suffix
            Sets a suffix appearing next to the number.
        valueformat
            Sets the value formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
        """

    def __init__(
        self, arg=None, font=None, prefix=None, suffix=None, valueformat=None, **kwargs
    ):
        """
        Construct a new Number object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.indicator.Number`
        font
            Set the font used to display main number
        prefix
            Sets a prefix appearing before the number.
        suffix
            Sets a suffix appearing next to the number.
        valueformat
            Sets the value formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.

        Returns
        -------
        Number
        """
        super().__init__("number")
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
The first argument to the plotly.graph_objs.indicator.Number
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.Number`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("font", arg, font)
        self._set_property("prefix", arg, prefix)
        self._set_property("suffix", arg, suffix)
        self._set_property("valueformat", arg, valueformat)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
