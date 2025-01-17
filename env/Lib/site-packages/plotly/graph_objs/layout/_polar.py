from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Polar(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.polar"
    _valid_props = {
        "angularaxis",
        "bargap",
        "barmode",
        "bgcolor",
        "domain",
        "gridshape",
        "hole",
        "radialaxis",
        "sector",
        "uirevision",
    }

    # angularaxis
    # -----------
    @property
    def angularaxis(self):
        """
        The 'angularaxis' property is an instance of AngularAxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.polar.AngularAxis`
          - A dict of string/value properties that will be passed
            to the AngularAxis constructor

            Supported dict properties:

                autotypenumbers
                    Using "strict" a numeric string in trace data
                    is not converted to a number. Using *convert
                    types* a numeric string in trace data may be
                    treated as a number during automatic axis
                    `type` detection. Defaults to
                    layout.autotypenumbers.
                categoryarray
                    Sets the order in which categories on this axis
                    appear. Only has an effect if `categoryorder`
                    is set to "array". Used with `categoryorder`.
                categoryarraysrc
                    Sets the source reference on Chart Studio Cloud
                    for `categoryarray`.
                categoryorder
                    Specifies the ordering logic for the case of
                    categorical variables. By default, plotly uses
                    "trace", which specifies the order that is
                    present in the data supplied. Set
                    `categoryorder` to *category ascending* or
                    *category descending* if order should be
                    determined by the alphanumerical order of the
                    category names. Set `categoryorder` to "array"
                    to derive the ordering from the attribute
                    `categoryarray`. If a category is not found in
                    the `categoryarray` array, the sorting behavior
                    for that attribute will be identical to the
                    "trace" mode. The unspecified categories will
                    follow the categories in `categoryarray`. Set
                    `categoryorder` to *total ascending* or *total
                    descending* if order should be determined by
                    the numerical order of the values. Similarly,
                    the order can be determined by the min, max,
                    sum, mean, geometric mean or median of all the
                    values.
                color
                    Sets default for all colors associated with
                    this axis all at once: line, font, tick, and
                    grid colors. Grid color is lightened by
                    blending this with the plot background
                    Individual pieces can override this.
                direction
                    Sets the direction corresponding to positive
                    angles.
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
                gridcolor
                    Sets the color of the grid lines.
                griddash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                gridwidth
                    Sets the width (in px) of the grid lines.
                hoverformat
                    Sets the hover text formatting rule using d3
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
                layer
                    Sets the layer on which this axis is displayed.
                    If *above traces*, this axis is displayed above
                    all the subplot's traces If *below traces*,
                    this axis is displayed below all the subplot's
                    traces, but above the grid lines. Useful when
                    used together with scatter-like traces with
                    `cliponaxis` set to False to show markers
                    and/or text nodes above this axis.
                linecolor
                    Sets the axis line color.
                linewidth
                    Sets the width (in px) of the axis line.
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
                period
                    Set the angular period. Has an effect only when
                    `angularaxis.type` is "category".
                rotation
                    Sets that start position (in degrees) of the
                    angular axis By default, polar subplots with
                    `direction` set to "counterclockwise" get a
                    `rotation` of 0 which corresponds to due East
                    (like what mathematicians prefer). In turn,
                    polar with `direction` set to "clockwise" get a
                    rotation of 90 which corresponds to due North
                    (like on a compass),
                separatethousands
                    If "true", even 4-digit integers are separated
                showexponent
                    If "all", all exponents are shown besides their
                    significands. If "first", only the exponent of
                    the first tick is shown. If "last", only the
                    exponent of the last tick is shown. If "none",
                    no exponents appear.
                showgrid
                    Determines whether or not grid lines are drawn.
                    If True, the grid lines are drawn at every tick
                    mark.
                showline
                    Determines whether or not a line bounding this
                    axis is drawn.
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
                thetaunit
                    Sets the format unit of the formatted "theta"
                    values. Has an effect only when
                    `angularaxis.type` is "linear".
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
                    Sets the tick font.
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
                    A tuple of :class:`plotly.graph_objects.layout.
                    polar.angularaxis.Tickformatstop` instances or
                    dicts with compatible properties
                tickformatstopdefaults
                    When used in a template (as layout.template.lay
                    out.polar.angularaxis.tickformatstopdefaults),
                    sets the default property values to use for
                    elements of
                    layout.polar.angularaxis.tickformatstops
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
                type
                    Sets the angular axis type. If "linear", set
                    `thetaunit` to determine the unit in which axis
                    value are shown. If *category, use `period` to
                    set the number of integer coordinates around
                    polar axis.
                uirevision
                    Controls persistence of user-driven changes in
                    axis `rotation`. Defaults to
                    `polar<N>.uirevision`.
                visible
                    A single toggle to hide the axis while
                    preserving interaction like dragging. Default
                    is true when a cheater plot is present on the
                    axis, otherwise false

        Returns
        -------
        plotly.graph_objs.layout.polar.AngularAxis
        """
        return self["angularaxis"]

    @angularaxis.setter
    def angularaxis(self, val):
        self["angularaxis"] = val

    # bargap
    # ------
    @property
    def bargap(self):
        """
        Sets the gap between bars of adjacent location coordinates.
        Values are unitless, they represent fractions of the minimum
        difference in bar positions in the data.

        The 'bargap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["bargap"]

    @bargap.setter
    def bargap(self, val):
        self["bargap"] = val

    # barmode
    # -------
    @property
    def barmode(self):
        """
        Determines how bars at the same location coordinate are
        displayed on the graph. With "stack", the bars are stacked on
        top of one another With "overlay", the bars are plotted over
        one another, you might need to reduce "opacity" to see multiple
        bars.

        The 'barmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['stack', 'overlay']

        Returns
        -------
        Any
        """
        return self["barmode"]

    @barmode.setter
    def barmode(self, val):
        self["barmode"] = val

    # bgcolor
    # -------
    @property
    def bgcolor(self):
        """
        Set the background color of the subplot

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

    # domain
    # ------
    @property
    def domain(self):
        """
        The 'domain' property is an instance of Domain
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.polar.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

            Supported dict properties:

                column
                    If there is a layout grid, use the domain for
                    this column in the grid for this polar subplot
                    .
                row
                    If there is a layout grid, use the domain for
                    this row in the grid for this polar subplot .
                x
                    Sets the horizontal domain of this polar
                    subplot (in plot fraction).
                y
                    Sets the vertical domain of this polar subplot
                    (in plot fraction).

        Returns
        -------
        plotly.graph_objs.layout.polar.Domain
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    # gridshape
    # ---------
    @property
    def gridshape(self):
        """
        Determines if the radial axis grid lines and angular axis line
        are drawn as "circular" sectors or as "linear" (polygon)
        sectors. Has an effect only when the angular axis has `type`
        "category". Note that `radialaxis.angle` is snapped to the
        angle of the closest vertex when `gridshape` is "circular" (so
        that radial axis scale is the same as the data scale).

        The 'gridshape' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['circular', 'linear']

        Returns
        -------
        Any
        """
        return self["gridshape"]

    @gridshape.setter
    def gridshape(self, val):
        self["gridshape"] = val

    # hole
    # ----
    @property
    def hole(self):
        """
        Sets the fraction of the radius to cut out of the polar
        subplot.

        The 'hole' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["hole"]

    @hole.setter
    def hole(self, val):
        self["hole"] = val

    # radialaxis
    # ----------
    @property
    def radialaxis(self):
        """
        The 'radialaxis' property is an instance of RadialAxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.polar.RadialAxis`
          - A dict of string/value properties that will be passed
            to the RadialAxis constructor

            Supported dict properties:

                angle
                    Sets the angle (in degrees) from which the
                    radial axis is drawn. Note that by default,
                    radial axis line on the theta=0 line
                    corresponds to a line pointing right (like what
                    mathematicians prefer). Defaults to the first
                    `polar.sector` angle.
                autorange
                    Determines whether or not the range of this
                    axis is computed in relation to the input data.
                    See `rangemode` for more info. If `range` is
                    provided and it has a value for both the lower
                    and upper bound, `autorange` is set to False.
                    Using "min" applies autorange only to set the
                    minimum. Using "max" applies autorange only to
                    set the maximum. Using *min reversed* applies
                    autorange only to set the minimum on a reversed
                    axis. Using *max reversed* applies autorange
                    only to set the maximum on a reversed axis.
                    Using "reversed" applies autorange on both ends
                    and reverses the axis direction.
                autorangeoptions
                    :class:`plotly.graph_objects.layout.polar.radia
                    laxis.Autorangeoptions` instance or dict with
                    compatible properties
                autotickangles
                    When `tickangle` is set to "auto", it will be
                    set to the first angle in this array that is
                    large enough to prevent label overlap.
                autotypenumbers
                    Using "strict" a numeric string in trace data
                    is not converted to a number. Using *convert
                    types* a numeric string in trace data may be
                    treated as a number during automatic axis
                    `type` detection. Defaults to
                    layout.autotypenumbers.
                calendar
                    Sets the calendar system to use for `range` and
                    `tick0` if this is a date axis. This does not
                    set the calendar for interpreting data on this
                    axis, that's specified in the trace or via the
                    global `layout.calendar`
                categoryarray
                    Sets the order in which categories on this axis
                    appear. Only has an effect if `categoryorder`
                    is set to "array". Used with `categoryorder`.
                categoryarraysrc
                    Sets the source reference on Chart Studio Cloud
                    for `categoryarray`.
                categoryorder
                    Specifies the ordering logic for the case of
                    categorical variables. By default, plotly uses
                    "trace", which specifies the order that is
                    present in the data supplied. Set
                    `categoryorder` to *category ascending* or
                    *category descending* if order should be
                    determined by the alphanumerical order of the
                    category names. Set `categoryorder` to "array"
                    to derive the ordering from the attribute
                    `categoryarray`. If a category is not found in
                    the `categoryarray` array, the sorting behavior
                    for that attribute will be identical to the
                    "trace" mode. The unspecified categories will
                    follow the categories in `categoryarray`. Set
                    `categoryorder` to *total ascending* or *total
                    descending* if order should be determined by
                    the numerical order of the values. Similarly,
                    the order can be determined by the min, max,
                    sum, mean, geometric mean or median of all the
                    values.
                color
                    Sets default for all colors associated with
                    this axis all at once: line, font, tick, and
                    grid colors. Grid color is lightened by
                    blending this with the plot background
                    Individual pieces can override this.
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
                gridcolor
                    Sets the color of the grid lines.
                griddash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                gridwidth
                    Sets the width (in px) of the grid lines.
                hoverformat
                    Sets the hover text formatting rule using d3
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
                layer
                    Sets the layer on which this axis is displayed.
                    If *above traces*, this axis is displayed above
                    all the subplot's traces If *below traces*,
                    this axis is displayed below all the subplot's
                    traces, but above the grid lines. Useful when
                    used together with scatter-like traces with
                    `cliponaxis` set to False to show markers
                    and/or text nodes above this axis.
                linecolor
                    Sets the axis line color.
                linewidth
                    Sets the width (in px) of the axis line.
                maxallowed
                    Determines the maximum range of this axis.
                minallowed
                    Determines the minimum range of this axis.
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
                    Sets the range of this axis. If the axis `type`
                    is "log", then you must take the log of your
                    desired range (e.g. to set the range from 1 to
                    100, set the range from 0 to 2). If the axis
                    `type` is "date", it should be date strings,
                    like date data, though Date objects and unix
                    milliseconds will be accepted and converted to
                    strings. If the axis `type` is "category", it
                    should be numbers, using the scale where each
                    category is assigned a serial number from zero
                    in the order it appears. Leaving either or both
                    elements `null` impacts the default
                    `autorange`.
                rangemode
                    If *tozero*`, the range extends to 0,
                    regardless of the input data If "nonnegative",
                    the range is non-negative, regardless of the
                    input data. If "normal", the range is computed
                    in relation to the extrema of the input data
                    (same behavior as for cartesian axes).
                separatethousands
                    If "true", even 4-digit integers are separated
                showexponent
                    If "all", all exponents are shown besides their
                    significands. If "first", only the exponent of
                    the first tick is shown. If "last", only the
                    exponent of the last tick is shown. If "none",
                    no exponents appear.
                showgrid
                    Determines whether or not grid lines are drawn.
                    If True, the grid lines are drawn at every tick
                    mark.
                showline
                    Determines whether or not a line bounding this
                    axis is drawn.
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
                side
                    Determines on which side of radial axis line
                    the tick and tick labels appear.
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
                    Sets the tick font.
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
                    A tuple of :class:`plotly.graph_objects.layout.
                    polar.radialaxis.Tickformatstop` instances or
                    dicts with compatible properties
                tickformatstopdefaults
                    When used in a template (as layout.template.lay
                    out.polar.radialaxis.tickformatstopdefaults),
                    sets the default property values to use for
                    elements of
                    layout.polar.radialaxis.tickformatstops
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
                title
                    :class:`plotly.graph_objects.layout.polar.radia
                    laxis.Title` instance or dict with compatible
                    properties
                titlefont
                    Deprecated: Please use
                    layout.polar.radialaxis.title.font instead.
                    Sets this axis' title font. Note that the
                    title's font used to be customized by the now
                    deprecated `titlefont` attribute.
                type
                    Sets the axis type. By default, plotly attempts
                    to determined the axis type by looking into the
                    data of the traces that referenced the axis in
                    question.
                uirevision
                    Controls persistence of user-driven changes in
                    axis `range`, `autorange`, `angle`, and `title`
                    if in `editable: true` configuration. Defaults
                    to `polar<N>.uirevision`.
                visible
                    A single toggle to hide the axis while
                    preserving interaction like dragging. Default
                    is true when a cheater plot is present on the
                    axis, otherwise false

        Returns
        -------
        plotly.graph_objs.layout.polar.RadialAxis
        """
        return self["radialaxis"]

    @radialaxis.setter
    def radialaxis(self, val):
        self["radialaxis"] = val

    # sector
    # ------
    @property
    def sector(self):
        """
            Sets angular span of this polar subplot with two angles (in
            degrees). Sector are assumed to be spanned in the
            counterclockwise direction with 0 corresponding to rightmost
            limit of the polar subplot.

            The 'sector' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'sector[0]' property is a number and may be specified as:
              - An int or float
        (1) The 'sector[1]' property is a number and may be specified as:
              - An int or float

            Returns
            -------
            list
        """
        return self["sector"]

    @sector.setter
    def sector(self, val):
        self["sector"] = val

    # uirevision
    # ----------
    @property
    def uirevision(self):
        """
        Controls persistence of user-driven changes in axis attributes,
        if not overridden in the individual axes. Defaults to
        `layout.uirevision`.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        angularaxis
            :class:`plotly.graph_objects.layout.polar.AngularAxis`
            instance or dict with compatible properties
        bargap
            Sets the gap between bars of adjacent location
            coordinates. Values are unitless, they represent
            fractions of the minimum difference in bar positions in
            the data.
        barmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "overlay", the bars
            are plotted over one another, you might need to reduce
            "opacity" to see multiple bars.
        bgcolor
            Set the background color of the subplot
        domain
            :class:`plotly.graph_objects.layout.polar.Domain`
            instance or dict with compatible properties
        gridshape
            Determines if the radial axis grid lines and angular
            axis line are drawn as "circular" sectors or as
            "linear" (polygon) sectors. Has an effect only when the
            angular axis has `type` "category". Note that
            `radialaxis.angle` is snapped to the angle of the
            closest vertex when `gridshape` is "circular" (so that
            radial axis scale is the same as the data scale).
        hole
            Sets the fraction of the radius to cut out of the polar
            subplot.
        radialaxis
            :class:`plotly.graph_objects.layout.polar.RadialAxis`
            instance or dict with compatible properties
        sector
            Sets angular span of this polar subplot with two angles
            (in degrees). Sector are assumed to be spanned in the
            counterclockwise direction with 0 corresponding to
            rightmost limit of the polar subplot.
        uirevision
            Controls persistence of user-driven changes in axis
            attributes, if not overridden in the individual axes.
            Defaults to `layout.uirevision`.
        """

    def __init__(
        self,
        arg=None,
        angularaxis=None,
        bargap=None,
        barmode=None,
        bgcolor=None,
        domain=None,
        gridshape=None,
        hole=None,
        radialaxis=None,
        sector=None,
        uirevision=None,
        **kwargs,
    ):
        """
        Construct a new Polar object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Polar`
        angularaxis
            :class:`plotly.graph_objects.layout.polar.AngularAxis`
            instance or dict with compatible properties
        bargap
            Sets the gap between bars of adjacent location
            coordinates. Values are unitless, they represent
            fractions of the minimum difference in bar positions in
            the data.
        barmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "overlay", the bars
            are plotted over one another, you might need to reduce
            "opacity" to see multiple bars.
        bgcolor
            Set the background color of the subplot
        domain
            :class:`plotly.graph_objects.layout.polar.Domain`
            instance or dict with compatible properties
        gridshape
            Determines if the radial axis grid lines and angular
            axis line are drawn as "circular" sectors or as
            "linear" (polygon) sectors. Has an effect only when the
            angular axis has `type` "category". Note that
            `radialaxis.angle` is snapped to the angle of the
            closest vertex when `gridshape` is "circular" (so that
            radial axis scale is the same as the data scale).
        hole
            Sets the fraction of the radius to cut out of the polar
            subplot.
        radialaxis
            :class:`plotly.graph_objects.layout.polar.RadialAxis`
            instance or dict with compatible properties
        sector
            Sets angular span of this polar subplot with two angles
            (in degrees). Sector are assumed to be spanned in the
            counterclockwise direction with 0 corresponding to
            rightmost limit of the polar subplot.
        uirevision
            Controls persistence of user-driven changes in axis
            attributes, if not overridden in the individual axes.
            Defaults to `layout.uirevision`.

        Returns
        -------
        Polar
        """
        super(Polar, self).__init__("polar")

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
The first argument to the plotly.graph_objs.layout.Polar
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Polar`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("angularaxis", None)
        _v = angularaxis if angularaxis is not None else _v
        if _v is not None:
            self["angularaxis"] = _v
        _v = arg.pop("bargap", None)
        _v = bargap if bargap is not None else _v
        if _v is not None:
            self["bargap"] = _v
        _v = arg.pop("barmode", None)
        _v = barmode if barmode is not None else _v
        if _v is not None:
            self["barmode"] = _v
        _v = arg.pop("bgcolor", None)
        _v = bgcolor if bgcolor is not None else _v
        if _v is not None:
            self["bgcolor"] = _v
        _v = arg.pop("domain", None)
        _v = domain if domain is not None else _v
        if _v is not None:
            self["domain"] = _v
        _v = arg.pop("gridshape", None)
        _v = gridshape if gridshape is not None else _v
        if _v is not None:
            self["gridshape"] = _v
        _v = arg.pop("hole", None)
        _v = hole if hole is not None else _v
        if _v is not None:
            self["hole"] = _v
        _v = arg.pop("radialaxis", None)
        _v = radialaxis if radialaxis is not None else _v
        if _v is not None:
            self["radialaxis"] = _v
        _v = arg.pop("sector", None)
        _v = sector if sector is not None else _v
        if _v is not None:
            self["sector"] = _v
        _v = arg.pop("uirevision", None)
        _v = uirevision if uirevision is not None else _v
        if _v is not None:
            self["uirevision"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
