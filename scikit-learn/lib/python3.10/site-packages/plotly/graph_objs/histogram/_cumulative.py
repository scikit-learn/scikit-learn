#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Cumulative(_BaseTraceHierarchyType):
    _parent_path_str = "histogram"
    _path_str = "histogram.cumulative"
    _valid_props = {"currentbin", "direction", "enabled"}

    @property
    def currentbin(self):
        """
        Only applies if cumulative is enabled. Sets whether the current
        bin is included, excluded, or has half of its value included in
        the current cumulative value. "include" is the default for
        compatibility with various other tools, however it introduces a
        half-bin bias to the results. "exclude" makes the opposite
        half-bin bias, and "half" removes it.

        The 'currentbin' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['include', 'exclude', 'half']

        Returns
        -------
        Any
        """
        return self["currentbin"]

    @currentbin.setter
    def currentbin(self, val):
        self["currentbin"] = val

    @property
    def direction(self):
        """
        Only applies if cumulative is enabled. If "increasing"
        (default) we sum all prior bins, so the result increases from
        left to right. If "decreasing" we sum later bins so the result
        decreases from left to right.

        The 'direction' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['increasing', 'decreasing']

        Returns
        -------
        Any
        """
        return self["direction"]

    @direction.setter
    def direction(self, val):
        self["direction"] = val

    @property
    def enabled(self):
        """
        If true, display the cumulative distribution by summing the
        binned values. Use the `direction` and `centralbin` attributes
        to tune the accumulation method. Note: in this mode, the
        "density" `histnorm` settings behave the same as their
        equivalents without "density": "" and "density" both rise to
        the number of data points, and "probability" and *probability
        density* both rise to the number of sample points.

        The 'enabled' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["enabled"]

    @enabled.setter
    def enabled(self, val):
        self["enabled"] = val

    @property
    def _prop_descriptions(self):
        return """\
        currentbin
            Only applies if cumulative is enabled. Sets whether the
            current bin is included, excluded, or has half of its
            value included in the current cumulative value.
            "include" is the default for compatibility with various
            other tools, however it introduces a half-bin bias to
            the results. "exclude" makes the opposite half-bin
            bias, and "half" removes it.
        direction
            Only applies if cumulative is enabled. If "increasing"
            (default) we sum all prior bins, so the result
            increases from left to right. If "decreasing" we sum
            later bins so the result decreases from left to right.
        enabled
            If true, display the cumulative distribution by summing
            the binned values. Use the `direction` and `centralbin`
            attributes to tune the accumulation method. Note: in
            this mode, the "density" `histnorm` settings behave the
            same as their equivalents without "density": "" and
            "density" both rise to the number of data points, and
            "probability" and *probability density* both rise to
            the number of sample points.
        """

    def __init__(
        self, arg=None, currentbin=None, direction=None, enabled=None, **kwargs
    ):
        """
        Construct a new Cumulative object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.histogram.Cumulative`
        currentbin
            Only applies if cumulative is enabled. Sets whether the
            current bin is included, excluded, or has half of its
            value included in the current cumulative value.
            "include" is the default for compatibility with various
            other tools, however it introduces a half-bin bias to
            the results. "exclude" makes the opposite half-bin
            bias, and "half" removes it.
        direction
            Only applies if cumulative is enabled. If "increasing"
            (default) we sum all prior bins, so the result
            increases from left to right. If "decreasing" we sum
            later bins so the result decreases from left to right.
        enabled
            If true, display the cumulative distribution by summing
            the binned values. Use the `direction` and `centralbin`
            attributes to tune the accumulation method. Note: in
            this mode, the "density" `histnorm` settings behave the
            same as their equivalents without "density": "" and
            "density" both rise to the number of data points, and
            "probability" and *probability density* both rise to
            the number of sample points.

        Returns
        -------
        Cumulative
        """
        super().__init__("cumulative")
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
The first argument to the plotly.graph_objs.histogram.Cumulative
constructor must be a dict or
an instance of :class:`plotly.graph_objs.histogram.Cumulative`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("currentbin", arg, currentbin)
        self._set_property("direction", arg, direction)
        self._set_property("enabled", arg, enabled)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
