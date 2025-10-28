#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class YAxis(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.xaxis.rangeslider"
    _path_str = "layout.xaxis.rangeslider.yaxis"
    _valid_props = {"range", "rangemode"}

    @property
    def range(self):
        """
            Sets the range of this axis for the rangeslider.

            The 'range' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'range[0]' property accepts values of any type
        (1) The 'range[1]' property accepts values of any type

            Returns
            -------
            list
        """
        return self["range"]

    @range.setter
    def range(self, val):
        self["range"] = val

    @property
    def rangemode(self):
        """
        Determines whether or not the range of this axis in the
        rangeslider use the same value than in the main plot when
        zooming in/out. If "auto", the autorange will be used. If
        "fixed", the `range` is used. If "match", the current range of
        the corresponding y-axis on the main subplot is used.

        The 'rangemode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['auto', 'fixed', 'match']

        Returns
        -------
        Any
        """
        return self["rangemode"]

    @rangemode.setter
    def rangemode(self, val):
        self["rangemode"] = val

    @property
    def _prop_descriptions(self):
        return """\
        range
            Sets the range of this axis for the rangeslider.
        rangemode
            Determines whether or not the range of this axis in the
            rangeslider use the same value than in the main plot
            when zooming in/out. If "auto", the autorange will be
            used. If "fixed", the `range` is used. If "match", the
            current range of the corresponding y-axis on the main
            subplot is used.
        """

    def __init__(self, arg=None, range=None, rangemode=None, **kwargs):
        """
        Construct a new YAxis object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.xaxis.r
            angeslider.YAxis`
        range
            Sets the range of this axis for the rangeslider.
        rangemode
            Determines whether or not the range of this axis in the
            rangeslider use the same value than in the main plot
            when zooming in/out. If "auto", the autorange will be
            used. If "fixed", the `range` is used. If "match", the
            current range of the corresponding y-axis on the main
            subplot is used.

        Returns
        -------
        YAxis
        """
        super().__init__("yaxis")
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
The first argument to the plotly.graph_objs.layout.xaxis.rangeslider.YAxis
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.xaxis.rangeslider.YAxis`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("range", arg, range)
        self._set_property("rangemode", arg, rangemode)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
