from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Gauge(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "indicator"
    _path_str = "indicator.gauge"
    _valid_props = {
        "axis",
        "bar",
        "bgcolor",
        "bordercolor",
        "borderwidth",
        "shape",
        "stepdefaults",
        "steps",
        "threshold",
    }

    # axis
    # ----
    @property
    def axis(self):
        """
        The 'axis' property is an instance of Axis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.gauge.Axis`
          - A dict of string/value properties that will be passed
            to the Axis constructor

            Supported dict properties:

                dtick
                    Sets the step in-between ticks on this axis.
                    Use with `tick0`. Must be a positive number, or
                    special strings available to "log" and "date"
                    axes. If the axis `type` is "log", then ticks
                    are set every 10^(n*dtick) where n is the tick
                    number. For example, to set a tick mark at 1,
                    10, 100, 1000, ... set dtick to 1. To set tick
                    marks at 1, 100, 10000, ... set dtick to 2. To
                    set tick marks at 1, 5, 25, 125, 625, 3125, ...
                    set dtick to log_10(5), or 0.69897000433. "log"
                    has several special values; "L<f>", where `f`
                    is a positive number, gives ticks linearly
                    spaced in value (but not position). For example
                    `tick0` = 0.1, `dtick` = "L0.5" will put ticks
                    at 0.1, 0.6, 1.1, 1.6 etc. To show powers of 10
                    plus small digits between, use "D1" (all
                    digits) or "D2" (only 2 and 5). `tick0` is
                    ignored for "D1" and "D2". If the axis `type`
                    is "date", then you must convert the time to
                    milliseconds. For example, to set the interval
                    between ticks to one day, set `dtick` to
                    86400000.0. "date" also has special values
                    "M<n>" gives ticks spaced by a number of
                    months. `n` must be a positive integer. To set
                    ticks on the 15th of every third month, set
                    `tick0` to "2000-01-15" and `dtick` to "M3". To
                    set ticks every 4 years, set `dtick` to "M48"
                exponentformat
                    Determines a formatting rule for the tick
                    exponents. For example, consider the number
                    1,000,000,000. If "none", it appears as
                    1,000,000,000. If "e", 1e+9. If "E", 1E+9. If
                    "power", 1x10^9 (with 9 in a super script). If
                    "SI", 1G. If "B", 1B.
                labelalias
                    Replacement text for specific tick or hover
                    labels. For example using {US: 'USA', CA:
                    'Canada'} changes US to USA and CA to Canada.
                    The labels we would have shown must match the
                    keys exactly, after adding any tickprefix or
                    ticksuffix. For negative numbers the minus sign
                    symbol used (U+2212) is wider than the regular
                    ascii dash. That means you need to use −1
                    instead of -1. labelalias can be used with any
                    axis type, and both keys (if needed) and values
                    (if desired) can include html-like tags or
                    MathJax.
                minexponent
                    Hide SI prefix for 10^n if |n| is below this
                    number. This only has an effect when
                    `tickformat` is "SI" or "B".
                nticks
                    Specifies the maximum number of ticks for the
                    particular axis. The actual number of ticks
                    will be chosen automatically to be less than or
                    equal to `nticks`. Has an effect only if
                    `tickmode` is set to "auto".
                range
                    Sets the range of this axis.
                separatethousands
                    If "true", even 4-digit integers are separated
                showexponent
                    If "all", all exponents are shown besides their
                    significands. If "first", only the exponent of
                    the first tick is shown. If "last", only the
                    exponent of the last tick is shown. If "none",
                    no exponents appear.
                showticklabels
                    Determines whether or not the tick labels are
                    drawn.
                showtickprefix
                    If "all", all tick labels are displayed with a
                    prefix. If "first", only the first tick is
                    displayed with a prefix. If "last", only the
                    last tick is displayed with a suffix. If
                    "none", tick prefixes are hidden.
                showticksuffix
                    Same as `showtickprefix` but for tick suffixes.
                tick0
                    Sets the placement of the first tick on this
                    axis. Use with `dtick`. If the axis `type` is
                    "log", then you must take the log of your
                    starting tick (e.g. to set the starting tick to
                    100, set the `tick0` to 2) except when
                    `dtick`=*L<f>* (see `dtick` for more info). If
                    the axis `type` is "date", it should be a date
                    string, like date data. If the axis `type` is
                    "category", it should be a number, using the
                    scale where each category is assigned a serial
                    number from zero in the order it appears.
                tickangle
                    Sets the angle of the tick labels with respect
                    to the horizontal. For example, a `tickangle`
                    of -90 draws the tick labels vertically.
                tickcolor
                    Sets the tick color.
                tickfont
                    Sets the color bar's tick label font
                tickformat
                    Sets the tick label formatting rule using d3
                    formatting mini-languages which are very
                    similar to those in Python. For numbers, see: h
                    ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                    format. And for dates see:
                    https://github.com/d3/d3-time-
                    format/tree/v2.2.3#locale_format. We add two
                    items to d3's date formatter: "%h" for half of
                    the year as a decimal number as well as "%{n}f"
                    for fractional seconds with n digits. For
                    example, *2016-10-13 09:15:23.456* with
                    tickformat "%H~%M~%S.%2f" would display
                    "09~15~23.46"
                tickformatstops
                    A tuple of :class:`plotly.graph_objects.indicat
                    or.gauge.axis.Tickformatstop` instances or
                    dicts with compatible properties
                tickformatstopdefaults
                    When used in a template (as layout.template.dat
                    a.indicator.gauge.axis.tickformatstopdefaults),
                    sets the default property values to use for
                    elements of
                    indicator.gauge.axis.tickformatstops
                ticklabelstep
                    Sets the spacing between tick labels as
                    compared to the spacing between ticks. A value
                    of 1 (default) means each tick gets a label. A
                    value of 2 means shows every 2nd label. A
                    larger value n means only every nth tick is
                    labeled. `tick0` determines which labels are
                    shown. Not implemented for axes with `type`
                    "log" or "multicategory", or when `tickmode` is
                    "array".
                ticklen
                    Sets the tick length (in px).
                tickmode
                    Sets the tick mode for this axis. If "auto",
                    the number of ticks is set via `nticks`. If
                    "linear", the placement of the ticks is
                    determined by a starting position `tick0` and a
                    tick step `dtick` ("linear" is the default
                    value if `tick0` and `dtick` are provided). If
                    "array", the placement of the ticks is set via
                    `tickvals` and the tick text is `ticktext`.
                    ("array" is the default value if `tickvals` is
                    provided).
                tickprefix
                    Sets a tick label prefix.
                ticks
                    Determines whether ticks are drawn or not. If
                    "", this axis' ticks are not drawn. If
                    "outside" ("inside"), this axis' are drawn
                    outside (inside) the axis lines.
                ticksuffix
                    Sets a tick label suffix.
                ticktext
                    Sets the text displayed at the ticks position
                    via `tickvals`. Only has an effect if
                    `tickmode` is set to "array". Used with
                    `tickvals`.
                ticktextsrc
                    Sets the source reference on Chart Studio Cloud
                    for `ticktext`.
                tickvals
                    Sets the values at which ticks on this axis
                    appear. Only has an effect if `tickmode` is set
                    to "array". Used with `ticktext`.
                tickvalssrc
                    Sets the source reference on Chart Studio Cloud
                    for `tickvals`.
                tickwidth
                    Sets the tick width (in px).
                visible
                    A single toggle to hide the axis while
                    preserving interaction like dragging. Default
                    is true when a cheater plot is present on the
                    axis, otherwise false

        Returns
        -------
        plotly.graph_objs.indicator.gauge.Axis
        """
        return self["axis"]

    @axis.setter
    def axis(self, val):
        self["axis"] = val

    # bar
    # ---
    @property
    def bar(self):
        """
        Set the appearance of the gauge's value

        The 'bar' property is an instance of Bar
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.gauge.Bar`
          - A dict of string/value properties that will be passed
            to the Bar constructor

            Supported dict properties:

                color
                    Sets the background color of the arc.
                line
                    :class:`plotly.graph_objects.indicator.gauge.ba
                    r.Line` instance or dict with compatible
                    properties
                thickness
                    Sets the thickness of the bar as a fraction of
                    the total thickness of the gauge.

        Returns
        -------
        plotly.graph_objs.indicator.gauge.Bar
        """
        return self["bar"]

    @bar.setter
    def bar(self, val):
        self["bar"] = val

    # bgcolor
    # -------
    @property
    def bgcolor(self):
        """
        Sets the gauge background color.

        The 'bgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color:
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, rebeccapurple, saddlebrown, salmon,
                sandybrown, seagreen, seashell, sienna, silver,
                skyblue, slateblue, slategray, slategrey, snow,
                springgreen, steelblue, tan, teal, thistle, tomato,
                turquoise, violet, wheat, white, whitesmoke,
                yellow, yellowgreen

        Returns
        -------
        str
        """
        return self["bgcolor"]

    @bgcolor.setter
    def bgcolor(self, val):
        self["bgcolor"] = val

    # bordercolor
    # -----------
    @property
    def bordercolor(self):
        """
        Sets the color of the border enclosing the gauge.

        The 'bordercolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color:
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, rebeccapurple, saddlebrown, salmon,
                sandybrown, seagreen, seashell, sienna, silver,
                skyblue, slateblue, slategray, slategrey, snow,
                springgreen, steelblue, tan, teal, thistle, tomato,
                turquoise, violet, wheat, white, whitesmoke,
                yellow, yellowgreen

        Returns
        -------
        str
        """
        return self["bordercolor"]

    @bordercolor.setter
    def bordercolor(self, val):
        self["bordercolor"] = val

    # borderwidth
    # -----------
    @property
    def borderwidth(self):
        """
        Sets the width (in px) of the border enclosing the gauge.

        The 'borderwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["borderwidth"]

    @borderwidth.setter
    def borderwidth(self, val):
        self["borderwidth"] = val

    # shape
    # -----
    @property
    def shape(self):
        """
        Set the shape of the gauge

        The 'shape' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['angular', 'bullet']

        Returns
        -------
        Any
        """
        return self["shape"]

    @shape.setter
    def shape(self, val):
        self["shape"] = val

    # steps
    # -----
    @property
    def steps(self):
        """
        The 'steps' property is a tuple of instances of
        Step that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.indicator.gauge.Step
          - A list or tuple of dicts of string/value properties that
            will be passed to the Step constructor

            Supported dict properties:

                color
                    Sets the background color of the arc.
                line
                    :class:`plotly.graph_objects.indicator.gauge.st
                    ep.Line` instance or dict with compatible
                    properties
                name
                    When used in a template, named items are
                    created in the output figure in addition to any
                    items the figure already has in this array. You
                    can modify these items in the output figure by
                    making your own item with `templateitemname`
                    matching this `name` alongside your
                    modifications (including `visible: false` or
                    `enabled: false` to hide it). Has no effect
                    outside of a template.
                range
                    Sets the range of this axis.
                templateitemname
                    Used to refer to a named item in this array in
                    the template. Named items from the template
                    will be created even without a matching item in
                    the input figure, but you can modify one by
                    making an item with `templateitemname` matching
                    its `name`, alongside your modifications
                    (including `visible: false` or `enabled: false`
                    to hide it). If there is no template or no
                    matching item, this item will be hidden unless
                    you explicitly show it with `visible: true`.
                thickness
                    Sets the thickness of the bar as a fraction of
                    the total thickness of the gauge.

        Returns
        -------
        tuple[plotly.graph_objs.indicator.gauge.Step]
        """
        return self["steps"]

    @steps.setter
    def steps(self, val):
        self["steps"] = val

    # stepdefaults
    # ------------
    @property
    def stepdefaults(self):
        """
        When used in a template (as
        layout.template.data.indicator.gauge.stepdefaults), sets the
        default property values to use for elements of
        indicator.gauge.steps

        The 'stepdefaults' property is an instance of Step
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.gauge.Step`
          - A dict of string/value properties that will be passed
            to the Step constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.indicator.gauge.Step
        """
        return self["stepdefaults"]

    @stepdefaults.setter
    def stepdefaults(self, val):
        self["stepdefaults"] = val

    # threshold
    # ---------
    @property
    def threshold(self):
        """
        The 'threshold' property is an instance of Threshold
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.gauge.Threshold`
          - A dict of string/value properties that will be passed
            to the Threshold constructor

            Supported dict properties:

                line
                    :class:`plotly.graph_objects.indicator.gauge.th
                    reshold.Line` instance or dict with compatible
                    properties
                thickness
                    Sets the thickness of the threshold line as a
                    fraction of the thickness of the gauge.
                value
                    Sets a treshold value drawn as a line.

        Returns
        -------
        plotly.graph_objs.indicator.gauge.Threshold
        """
        return self["threshold"]

    @threshold.setter
    def threshold(self, val):
        self["threshold"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        axis
            :class:`plotly.graph_objects.indicator.gauge.Axis`
            instance or dict with compatible properties
        bar
            Set the appearance of the gauge's value
        bgcolor
            Sets the gauge background color.
        bordercolor
            Sets the color of the border enclosing the gauge.
        borderwidth
            Sets the width (in px) of the border enclosing the
            gauge.
        shape
            Set the shape of the gauge
        steps
            A tuple of
            :class:`plotly.graph_objects.indicator.gauge.Step`
            instances or dicts with compatible properties
        stepdefaults
            When used in a template (as
            layout.template.data.indicator.gauge.stepdefaults),
            sets the default property values to use for elements of
            indicator.gauge.steps
        threshold
            :class:`plotly.graph_objects.indicator.gauge.Threshold`
            instance or dict with compatible properties
        """

    def __init__(
        self,
        arg=None,
        axis=None,
        bar=None,
        bgcolor=None,
        bordercolor=None,
        borderwidth=None,
        shape=None,
        steps=None,
        stepdefaults=None,
        threshold=None,
        **kwargs,
    ):
        """
        Construct a new Gauge object

        The gauge of the Indicator plot.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.indicator.Gauge`
        axis
            :class:`plotly.graph_objects.indicator.gauge.Axis`
            instance or dict with compatible properties
        bar
            Set the appearance of the gauge's value
        bgcolor
            Sets the gauge background color.
        bordercolor
            Sets the color of the border enclosing the gauge.
        borderwidth
            Sets the width (in px) of the border enclosing the
            gauge.
        shape
            Set the shape of the gauge
        steps
            A tuple of
            :class:`plotly.graph_objects.indicator.gauge.Step`
            instances or dicts with compatible properties
        stepdefaults
            When used in a template (as
            layout.template.data.indicator.gauge.stepdefaults),
            sets the default property values to use for elements of
            indicator.gauge.steps
        threshold
            :class:`plotly.graph_objects.indicator.gauge.Threshold`
            instance or dict with compatible properties

        Returns
        -------
        Gauge
        """
        super(Gauge, self).__init__("gauge")

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
The first argument to the plotly.graph_objs.indicator.Gauge
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.Gauge`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("axis", None)
        _v = axis if axis is not None else _v
        if _v is not None:
            self["axis"] = _v
        _v = arg.pop("bar", None)
        _v = bar if bar is not None else _v
        if _v is not None:
            self["bar"] = _v
        _v = arg.pop("bgcolor", None)
        _v = bgcolor if bgcolor is not None else _v
        if _v is not None:
            self["bgcolor"] = _v
        _v = arg.pop("bordercolor", None)
        _v = bordercolor if bordercolor is not None else _v
        if _v is not None:
            self["bordercolor"] = _v
        _v = arg.pop("borderwidth", None)
        _v = borderwidth if borderwidth is not None else _v
        if _v is not None:
            self["borderwidth"] = _v
        _v = arg.pop("shape", None)
        _v = shape if shape is not None else _v
        if _v is not None:
            self["shape"] = _v
        _v = arg.pop("steps", None)
        _v = steps if steps is not None else _v
        if _v is not None:
            self["steps"] = _v
        _v = arg.pop("stepdefaults", None)
        _v = stepdefaults if stepdefaults is not None else _v
        if _v is not None:
            self["stepdefaults"] = _v
        _v = arg.pop("threshold", None)
        _v = threshold if threshold is not None else _v
        if _v is not None:
            self["threshold"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
