from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Delta(_BaseTraceHierarchyType):

    # class properties
    # --------------------
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

    # decreasing
    # ----------
    @property
    def decreasing(self):
        """
        The 'decreasing' property is an instance of Decreasing
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.delta.Decreasing`
          - A dict of string/value properties that will be passed
            to the Decreasing constructor

            Supported dict properties:

                color
                    Sets the color for increasing value.
                symbol
                    Sets the symbol to display for increasing value

        Returns
        -------
        plotly.graph_objs.indicator.delta.Decreasing
        """
        return self["decreasing"]

    @decreasing.setter
    def decreasing(self, val):
        self["decreasing"] = val

    # font
    # ----
    @property
    def font(self):
        """
        Set the font used to display the delta

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.delta.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

            Supported dict properties:

                color

                family
                    HTML font family - the typeface that will be
                    applied by the web browser. The web browser
                    will only be able to apply a font if it is
                    available on the system which it operates.
                    Provide multiple font families, separated by
                    commas, to indicate the preference in which to
                    apply fonts if they aren't available on the
                    system. The Chart Studio Cloud (at
                    https://chart-studio.plotly.com or on-premise)
                    generates images on a server, where only a
                    select number of fonts are installed and
                    supported. These include "Arial", "Balto",
                    "Courier New", "Droid Sans", "Droid Serif",
                    "Droid Sans Mono", "Gravitas One", "Old
                    Standard TT", "Open Sans", "Overpass", "PT Sans
                    Narrow", "Raleway", "Times New Roman".
                lineposition
                    Sets the kind of decoration line(s) with text,
                    such as an "under", "over" or "through" as well
                    as combinations e.g. "under+over", etc.
                shadow
                    Sets the shape and color of the shadow behind
                    text. "auto" places minimal shadow and applies
                    contrast text font color. See
                    https://developer.mozilla.org/en-
                    US/docs/Web/CSS/text-shadow for additional
                    options.
                size

                style
                    Sets whether a font should be styled with a
                    normal or italic face from its family.
                textcase
                    Sets capitalization of text. It can be used to
                    make text appear in all-uppercase or all-
                    lowercase, or with each word capitalized.
                variant
                    Sets the variant of the font.
                weight
                    Sets the weight (or boldness) of the font.

        Returns
        -------
        plotly.graph_objs.indicator.delta.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # increasing
    # ----------
    @property
    def increasing(self):
        """
        The 'increasing' property is an instance of Increasing
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.delta.Increasing`
          - A dict of string/value properties that will be passed
            to the Increasing constructor

            Supported dict properties:

                color
                    Sets the color for increasing value.
                symbol
                    Sets the symbol to display for increasing value

        Returns
        -------
        plotly.graph_objs.indicator.delta.Increasing
        """
        return self["increasing"]

    @increasing.setter
    def increasing(self, val):
        self["increasing"] = val

    # position
    # --------
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

    # prefix
    # ------
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

    # reference
    # ---------
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

    # relative
    # --------
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

    # suffix
    # ------
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

    # valueformat
    # -----------
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

    # Self properties description
    # ---------------------------
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
        super(Delta, self).__init__("delta")

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
The first argument to the plotly.graph_objs.indicator.Delta
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.Delta`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("decreasing", None)
        _v = decreasing if decreasing is not None else _v
        if _v is not None:
            self["decreasing"] = _v
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("increasing", None)
        _v = increasing if increasing is not None else _v
        if _v is not None:
            self["increasing"] = _v
        _v = arg.pop("position", None)
        _v = position if position is not None else _v
        if _v is not None:
            self["position"] = _v
        _v = arg.pop("prefix", None)
        _v = prefix if prefix is not None else _v
        if _v is not None:
            self["prefix"] = _v
        _v = arg.pop("reference", None)
        _v = reference if reference is not None else _v
        if _v is not None:
            self["reference"] = _v
        _v = arg.pop("relative", None)
        _v = relative if relative is not None else _v
        if _v is not None:
            self["relative"] = _v
        _v = arg.pop("suffix", None)
        _v = suffix if suffix is not None else _v
        if _v is not None:
            self["suffix"] = _v
        _v = arg.pop("valueformat", None)
        _v = valueformat if valueformat is not None else _v
        if _v is not None:
            self["valueformat"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
