from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "scatter"
    _path_str = "scatter.marker"
    _valid_props = {
        "angle",
        "angleref",
        "anglesrc",
        "autocolorscale",
        "cauto",
        "cmax",
        "cmid",
        "cmin",
        "color",
        "coloraxis",
        "colorbar",
        "colorscale",
        "colorsrc",
        "gradient",
        "line",
        "maxdisplayed",
        "opacity",
        "opacitysrc",
        "reversescale",
        "showscale",
        "size",
        "sizemin",
        "sizemode",
        "sizeref",
        "sizesrc",
        "standoff",
        "standoffsrc",
        "symbol",
        "symbolsrc",
    }

    # angle
    # -----
    @property
    def angle(self):
        """
        Sets the marker angle in respect to `angleref`.

        The 'angle' property is a angle (in degrees) that may be
        specified as a number between -180 and 180, or a list, numpy array or other iterable thereof.
        Numeric values outside this range are converted to the equivalent value
        (e.g. 270 is converted to -90).

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["angle"]

    @angle.setter
    def angle(self, val):
        self["angle"] = val

    # angleref
    # --------
    @property
    def angleref(self):
        """
        Sets the reference for marker angle. With "previous", angle 0
        points along the line from the previous point to this one. With
        "up", angle 0 points toward the top of the screen.

        The 'angleref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['previous', 'up']

        Returns
        -------
        Any
        """
        return self["angleref"]

    @angleref.setter
    def angleref(self, val):
        self["angleref"] = val

    # anglesrc
    # --------
    @property
    def anglesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `angle`.

        The 'anglesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["anglesrc"]

    @anglesrc.setter
    def anglesrc(self, val):
        self["anglesrc"] = val

    # autocolorscale
    # --------------
    @property
    def autocolorscale(self):
        """
        Determines whether the colorscale is a default palette
        (`autocolorscale: true`) or the palette determined by
        `marker.colorscale`. Has an effect only if in `marker.color` is
        set to a numerical array. In case `colorscale` is unspecified
        or `autocolorscale` is true, the default palette will be chosen
        according to whether numbers in the `color` array are all
        positive, all negative or mixed.

        The 'autocolorscale' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["autocolorscale"]

    @autocolorscale.setter
    def autocolorscale(self, val):
        self["autocolorscale"] = val

    # cauto
    # -----
    @property
    def cauto(self):
        """
        Determines whether or not the color domain is computed with
        respect to the input data (here in `marker.color`) or the
        bounds set in `marker.cmin` and `marker.cmax` Has an effect
        only if in `marker.color` is set to a numerical array. Defaults
        to `false` when `marker.cmin` and `marker.cmax` are set by the
        user.

        The 'cauto' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["cauto"]

    @cauto.setter
    def cauto(self, val):
        self["cauto"] = val

    # cmax
    # ----
    @property
    def cmax(self):
        """
        Sets the upper bound of the color domain. Has an effect only if
        in `marker.color` is set to a numerical array. Value should
        have the same units as in `marker.color` and if set,
        `marker.cmin` must be set as well.

        The 'cmax' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["cmax"]

    @cmax.setter
    def cmax(self, val):
        self["cmax"] = val

    # cmid
    # ----
    @property
    def cmid(self):
        """
        Sets the mid-point of the color domain by scaling `marker.cmin`
        and/or `marker.cmax` to be equidistant to this point. Has an
        effect only if in `marker.color` is set to a numerical array.
        Value should have the same units as in `marker.color`. Has no
        effect when `marker.cauto` is `false`.

        The 'cmid' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["cmid"]

    @cmid.setter
    def cmid(self, val):
        self["cmid"] = val

    # cmin
    # ----
    @property
    def cmin(self):
        """
        Sets the lower bound of the color domain. Has an effect only if
        in `marker.color` is set to a numerical array. Value should
        have the same units as in `marker.color` and if set,
        `marker.cmax` must be set as well.

        The 'cmin' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["cmin"]

    @cmin.setter
    def cmin(self, val):
        self["cmin"] = val

    # color
    # -----
    @property
    def color(self):
        """
        Sets the marker color. It accepts either a specific color or an
        array of numbers that are mapped to the colorscale relative to
        the max and min values of the array or relative to
        `marker.cmin` and `marker.cmax` if set.

        The 'color' property is a color and may be specified as:
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
          - A number that will be interpreted as a color
            according to scatter.marker.colorscale
          - A list or array of any of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    # coloraxis
    # ---------
    @property
    def coloraxis(self):
        """
        Sets a reference to a shared color axis. References to these
        shared color axes are "coloraxis", "coloraxis2", "coloraxis3",
        etc. Settings for these shared color axes are set in the
        layout, under `layout.coloraxis`, `layout.coloraxis2`, etc.
        Note that multiple color scales can be linked to the same color
        axis.

        The 'coloraxis' property is an identifier of a particular
        subplot, of type 'coloraxis', that may be specified as the string 'coloraxis'
        optionally followed by an integer >= 1
        (e.g. 'coloraxis', 'coloraxis1', 'coloraxis2', 'coloraxis3', etc.)

        Returns
        -------
        str
        """
        return self["coloraxis"]

    @coloraxis.setter
    def coloraxis(self, val):
        self["coloraxis"] = val

    # colorbar
    # --------
    @property
    def colorbar(self):
        """
        The 'colorbar' property is an instance of ColorBar
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatter.marker.ColorBar`
          - A dict of string/value properties that will be passed
            to the ColorBar constructor

            Supported dict properties:

                bgcolor
                    Sets the color of padded area.
                bordercolor
                    Sets the axis line color.
                borderwidth
                    Sets the width (in px) or the border enclosing
                    this color bar.
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
                len
                    Sets the length of the color bar This measure
                    excludes the padding of both ends. That is, the
                    color bar length is this length minus the
                    padding on both ends.
                lenmode
                    Determines whether this color bar's length
                    (i.e. the measure in the color variation
                    direction) is set in units of plot "fraction"
                    or in *pixels. Use `len` to set the value.
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
                orientation
                    Sets the orientation of the colorbar.
                outlinecolor
                    Sets the axis line color.
                outlinewidth
                    Sets the width (in px) of the axis line.
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
                thickness
                    Sets the thickness of the color bar This
                    measure excludes the size of the padding, ticks
                    and labels.
                thicknessmode
                    Determines whether this color bar's thickness
                    (i.e. the measure in the constant color
                    direction) is set in units of plot "fraction"
                    or in "pixels". Use `thickness` to set the
                    value.
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
                    A tuple of :class:`plotly.graph_objects.scatter
                    .marker.colorbar.Tickformatstop` instances or
                    dicts with compatible properties
                tickformatstopdefaults
                    When used in a template (as layout.template.dat
                    a.scatter.marker.colorbar.tickformatstopdefault
                    s), sets the default property values to use for
                    elements of
                    scatter.marker.colorbar.tickformatstops
                ticklabeloverflow
                    Determines how we handle tick labels that would
                    overflow either the graph div or the domain of
                    the axis. The default value for inside tick
                    labels is *hide past domain*. In other cases
                    the default is *hide past div*.
                ticklabelposition
                    Determines where tick labels are drawn relative
                    to the ticks. Left and right options are used
                    when `orientation` is "h", top and bottom when
                    `orientation` is "v".
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
                    :class:`plotly.graph_objects.scatter.marker.col
                    orbar.Title` instance or dict with compatible
                    properties
                titlefont
                    Deprecated: Please use
                    scatter.marker.colorbar.title.font instead.
                    Sets this color bar's title font. Note that the
                    title's font used to be set by the now
                    deprecated `titlefont` attribute.
                titleside
                    Deprecated: Please use
                    scatter.marker.colorbar.title.side instead.
                    Determines the location of color bar's title
                    with respect to the color bar. Defaults to
                    "top" when `orientation` if "v" and  defaults
                    to "right" when `orientation` if "h". Note that
                    the title's location used to be set by the now
                    deprecated `titleside` attribute.
                x
                    Sets the x position with respect to `xref` of
                    the color bar (in plot fraction). When `xref`
                    is "paper", defaults to 1.02 when `orientation`
                    is "v" and 0.5 when `orientation` is "h". When
                    `xref` is "container", defaults to 1 when
                    `orientation` is "v" and 0.5 when `orientation`
                    is "h". Must be between 0 and 1 if `xref` is
                    "container" and between "-2" and 3 if `xref` is
                    "paper".
                xanchor
                    Sets this color bar's horizontal position
                    anchor. This anchor binds the `x` position to
                    the "left", "center" or "right" of the color
                    bar. Defaults to "left" when `orientation` is
                    "v" and "center" when `orientation` is "h".
                xpad
                    Sets the amount of padding (in px) along the x
                    direction.
                xref
                    Sets the container `x` refers to. "container"
                    spans the entire `width` of the plot. "paper"
                    refers to the width of the plotting area only.
                y
                    Sets the y position with respect to `yref` of
                    the color bar (in plot fraction). When `yref`
                    is "paper", defaults to 0.5 when `orientation`
                    is "v" and 1.02 when `orientation` is "h". When
                    `yref` is "container", defaults to 0.5 when
                    `orientation` is "v" and 1 when `orientation`
                    is "h". Must be between 0 and 1 if `yref` is
                    "container" and between "-2" and 3 if `yref` is
                    "paper".
                yanchor
                    Sets this color bar's vertical position anchor
                    This anchor binds the `y` position to the
                    "top", "middle" or "bottom" of the color bar.
                    Defaults to "middle" when `orientation` is "v"
                    and "bottom" when `orientation` is "h".
                ypad
                    Sets the amount of padding (in px) along the y
                    direction.
                yref
                    Sets the container `y` refers to. "container"
                    spans the entire `height` of the plot. "paper"
                    refers to the height of the plotting area only.

        Returns
        -------
        plotly.graph_objs.scatter.marker.ColorBar
        """
        return self["colorbar"]

    @colorbar.setter
    def colorbar(self, val):
        self["colorbar"] = val

    # colorscale
    # ----------
    @property
    def colorscale(self):
        """
        Sets the colorscale. Has an effect only if in `marker.color` is
        set to a numerical array. The colorscale must be an array
        containing arrays mapping a normalized value to an rgb, rgba,
        hex, hsl, hsv, or named color string. At minimum, a mapping for
        the lowest (0) and highest (1) values are required. For
        example, `[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`. To
        control the bounds of the colorscale in color space, use
        `marker.cmin` and `marker.cmax`. Alternatively, `colorscale`
        may be a palette name string of the following list: Blackbody,B
        luered,Blues,Cividis,Earth,Electric,Greens,Greys,Hot,Jet,Picnic
        ,Portland,Rainbow,RdBu,Reds,Viridis,YlGnBu,YlOrRd.

        The 'colorscale' property is a colorscale and may be
        specified as:
          - A list of colors that will be spaced evenly to create the colorscale.
            Many predefined colorscale lists are included in the sequential, diverging,
            and cyclical modules in the plotly.colors package.
          - A list of 2-element lists where the first element is the
            normalized color level value (starting at 0 and ending at 1),
            and the second item is a valid color string.
            (e.g. [[0, 'green'], [0.5, 'red'], [1.0, 'rgb(0, 0, 255)']])
          - One of the following named colorscales:
                ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
                 'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
                 'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
                 'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
                 'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
                 'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
                 'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
                 'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
                 'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
                 'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
                 'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
                 'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
                 'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
                 'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr',
                 'ylorrd'].
            Appending '_r' to a named colorscale reverses it.

        Returns
        -------
        str
        """
        return self["colorscale"]

    @colorscale.setter
    def colorscale(self, val):
        self["colorscale"] = val

    # colorsrc
    # --------
    @property
    def colorsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `color`.

        The 'colorsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["colorsrc"]

    @colorsrc.setter
    def colorsrc(self, val):
        self["colorsrc"] = val

    # gradient
    # --------
    @property
    def gradient(self):
        """
        The 'gradient' property is an instance of Gradient
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatter.marker.Gradient`
          - A dict of string/value properties that will be passed
            to the Gradient constructor

            Supported dict properties:

                color
                    Sets the final color of the gradient fill: the
                    center color for radial, the right for
                    horizontal, or the bottom for vertical.
                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.
                type
                    Sets the type of gradient used to fill the
                    markers
                typesrc
                    Sets the source reference on Chart Studio Cloud
                    for `type`.

        Returns
        -------
        plotly.graph_objs.scatter.marker.Gradient
        """
        return self["gradient"]

    @gradient.setter
    def gradient(self, val):
        self["gradient"] = val

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatter.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                autocolorscale
                    Determines whether the colorscale is a default
                    palette (`autocolorscale: true`) or the palette
                    determined by `marker.line.colorscale`. Has an
                    effect only if in `marker.line.color` is set to
                    a numerical array. In case `colorscale` is
                    unspecified or `autocolorscale` is true, the
                    default palette will be chosen according to
                    whether numbers in the `color` array are all
                    positive, all negative or mixed.
                cauto
                    Determines whether or not the color domain is
                    computed with respect to the input data (here
                    in `marker.line.color`) or the bounds set in
                    `marker.line.cmin` and `marker.line.cmax` Has
                    an effect only if in `marker.line.color` is set
                    to a numerical array. Defaults to `false` when
                    `marker.line.cmin` and `marker.line.cmax` are
                    set by the user.
                cmax
                    Sets the upper bound of the color domain. Has
                    an effect only if in `marker.line.color` is set
                    to a numerical array. Value should have the
                    same units as in `marker.line.color` and if
                    set, `marker.line.cmin` must be set as well.
                cmid
                    Sets the mid-point of the color domain by
                    scaling `marker.line.cmin` and/or
                    `marker.line.cmax` to be equidistant to this
                    point. Has an effect only if in
                    `marker.line.color` is set to a numerical
                    array. Value should have the same units as in
                    `marker.line.color`. Has no effect when
                    `marker.line.cauto` is `false`.
                cmin
                    Sets the lower bound of the color domain. Has
                    an effect only if in `marker.line.color` is set
                    to a numerical array. Value should have the
                    same units as in `marker.line.color` and if
                    set, `marker.line.cmax` must be set as well.
                color
                    Sets the marker.line color. It accepts either a
                    specific color or an array of numbers that are
                    mapped to the colorscale relative to the max
                    and min values of the array or relative to
                    `marker.line.cmin` and `marker.line.cmax` if
                    set.
                coloraxis
                    Sets a reference to a shared color axis.
                    References to these shared color axes are
                    "coloraxis", "coloraxis2", "coloraxis3", etc.
                    Settings for these shared color axes are set in
                    the layout, under `layout.coloraxis`,
                    `layout.coloraxis2`, etc. Note that multiple
                    color scales can be linked to the same color
                    axis.
                colorscale
                    Sets the colorscale. Has an effect only if in
                    `marker.line.color` is set to a numerical
                    array. The colorscale must be an array
                    containing arrays mapping a normalized value to
                    an rgb, rgba, hex, hsl, hsv, or named color
                    string. At minimum, a mapping for the lowest
                    (0) and highest (1) values are required. For
                    example, `[[0, 'rgb(0,0,255)'], [1,
                    'rgb(255,0,0)']]`. To control the bounds of the
                    colorscale in color space, use
                    `marker.line.cmin` and `marker.line.cmax`.
                    Alternatively, `colorscale` may be a palette
                    name string of the following list: Blackbody,Bl
                    uered,Blues,Cividis,Earth,Electric,Greens,Greys
                    ,Hot,Jet,Picnic,Portland,Rainbow,RdBu,Reds,Viri
                    dis,YlGnBu,YlOrRd.
                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.
                reversescale
                    Reverses the color mapping if true. Has an
                    effect only if in `marker.line.color` is set to
                    a numerical array. If true, `marker.line.cmin`
                    will correspond to the last color in the array
                    and `marker.line.cmax` will correspond to the
                    first color.
                width
                    Sets the width (in px) of the lines bounding
                    the marker points.
                widthsrc
                    Sets the source reference on Chart Studio Cloud
                    for `width`.

        Returns
        -------
        plotly.graph_objs.scatter.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # maxdisplayed
    # ------------
    @property
    def maxdisplayed(self):
        """
        Sets a maximum number of points to be drawn on the graph. 0
        corresponds to no limit.

        The 'maxdisplayed' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["maxdisplayed"]

    @maxdisplayed.setter
    def maxdisplayed(self, val):
        self["maxdisplayed"] = val

    # opacity
    # -------
    @property
    def opacity(self):
        """
        Sets the marker opacity.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    # opacitysrc
    # ----------
    @property
    def opacitysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `opacity`.

        The 'opacitysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["opacitysrc"]

    @opacitysrc.setter
    def opacitysrc(self, val):
        self["opacitysrc"] = val

    # reversescale
    # ------------
    @property
    def reversescale(self):
        """
        Reverses the color mapping if true. Has an effect only if in
        `marker.color` is set to a numerical array. If true,
        `marker.cmin` will correspond to the last color in the array
        and `marker.cmax` will correspond to the first color.

        The 'reversescale' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["reversescale"]

    @reversescale.setter
    def reversescale(self, val):
        self["reversescale"] = val

    # showscale
    # ---------
    @property
    def showscale(self):
        """
        Determines whether or not a colorbar is displayed for this
        trace. Has an effect only if in `marker.color` is set to a
        numerical array.

        The 'showscale' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showscale"]

    @showscale.setter
    def showscale(self, val):
        self["showscale"] = val

    # size
    # ----
    @property
    def size(self):
        """
        Sets the marker size (in px).

        The 'size' property is a number and may be specified as:
          - An int or float in the interval [0, inf]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    # sizemin
    # -------
    @property
    def sizemin(self):
        """
        Has an effect only if `marker.size` is set to a numerical
        array. Sets the minimum size (in px) of the rendered marker
        points.

        The 'sizemin' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["sizemin"]

    @sizemin.setter
    def sizemin(self, val):
        self["sizemin"] = val

    # sizemode
    # --------
    @property
    def sizemode(self):
        """
        Has an effect only if `marker.size` is set to a numerical
        array. Sets the rule for which the data in `size` is converted
        to pixels.

        The 'sizemode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['diameter', 'area']

        Returns
        -------
        Any
        """
        return self["sizemode"]

    @sizemode.setter
    def sizemode(self, val):
        self["sizemode"] = val

    # sizeref
    # -------
    @property
    def sizeref(self):
        """
        Has an effect only if `marker.size` is set to a numerical
        array. Sets the scale factor used to determine the rendered
        size of marker points. Use with `sizemin` and `sizemode`.

        The 'sizeref' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["sizeref"]

    @sizeref.setter
    def sizeref(self, val):
        self["sizeref"] = val

    # sizesrc
    # -------
    @property
    def sizesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `size`.

        The 'sizesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["sizesrc"]

    @sizesrc.setter
    def sizesrc(self, val):
        self["sizesrc"] = val

    # standoff
    # --------
    @property
    def standoff(self):
        """
        Moves the marker away from the data point in the direction of
        `angle` (in px). This can be useful for example if you have
        another marker at this location and you want to point an
        arrowhead marker at it.

        The 'standoff' property is a number and may be specified as:
          - An int or float in the interval [0, inf]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["standoff"]

    @standoff.setter
    def standoff(self, val):
        self["standoff"] = val

    # standoffsrc
    # -----------
    @property
    def standoffsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `standoff`.

        The 'standoffsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["standoffsrc"]

    @standoffsrc.setter
    def standoffsrc(self, val):
        self["standoffsrc"] = val

    # symbol
    # ------
    @property
    def symbol(self):
        """
        Sets the marker symbol type. Adding 100 is equivalent to
        appending "-open" to a symbol name. Adding 200 is equivalent to
        appending "-dot" to a symbol name. Adding 300 is equivalent to
        appending "-open-dot" or "dot-open" to a symbol name.

        The 'symbol' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [0, '0', 'circle', 100, '100', 'circle-open', 200, '200',
                'circle-dot', 300, '300', 'circle-open-dot', 1, '1',
                'square', 101, '101', 'square-open', 201, '201',
                'square-dot', 301, '301', 'square-open-dot', 2, '2',
                'diamond', 102, '102', 'diamond-open', 202, '202',
                'diamond-dot', 302, '302', 'diamond-open-dot', 3, '3',
                'cross', 103, '103', 'cross-open', 203, '203',
                'cross-dot', 303, '303', 'cross-open-dot', 4, '4', 'x',
                104, '104', 'x-open', 204, '204', 'x-dot', 304, '304',
                'x-open-dot', 5, '5', 'triangle-up', 105, '105',
                'triangle-up-open', 205, '205', 'triangle-up-dot', 305,
                '305', 'triangle-up-open-dot', 6, '6', 'triangle-down',
                106, '106', 'triangle-down-open', 206, '206',
                'triangle-down-dot', 306, '306', 'triangle-down-open-dot',
                7, '7', 'triangle-left', 107, '107', 'triangle-left-open',
                207, '207', 'triangle-left-dot', 307, '307',
                'triangle-left-open-dot', 8, '8', 'triangle-right', 108,
                '108', 'triangle-right-open', 208, '208',
                'triangle-right-dot', 308, '308',
                'triangle-right-open-dot', 9, '9', 'triangle-ne', 109,
                '109', 'triangle-ne-open', 209, '209', 'triangle-ne-dot',
                309, '309', 'triangle-ne-open-dot', 10, '10',
                'triangle-se', 110, '110', 'triangle-se-open', 210, '210',
                'triangle-se-dot', 310, '310', 'triangle-se-open-dot', 11,
                '11', 'triangle-sw', 111, '111', 'triangle-sw-open', 211,
                '211', 'triangle-sw-dot', 311, '311',
                'triangle-sw-open-dot', 12, '12', 'triangle-nw', 112,
                '112', 'triangle-nw-open', 212, '212', 'triangle-nw-dot',
                312, '312', 'triangle-nw-open-dot', 13, '13', 'pentagon',
                113, '113', 'pentagon-open', 213, '213', 'pentagon-dot',
                313, '313', 'pentagon-open-dot', 14, '14', 'hexagon', 114,
                '114', 'hexagon-open', 214, '214', 'hexagon-dot', 314,
                '314', 'hexagon-open-dot', 15, '15', 'hexagon2', 115,
                '115', 'hexagon2-open', 215, '215', 'hexagon2-dot', 315,
                '315', 'hexagon2-open-dot', 16, '16', 'octagon', 116,
                '116', 'octagon-open', 216, '216', 'octagon-dot', 316,
                '316', 'octagon-open-dot', 17, '17', 'star', 117, '117',
                'star-open', 217, '217', 'star-dot', 317, '317',
                'star-open-dot', 18, '18', 'hexagram', 118, '118',
                'hexagram-open', 218, '218', 'hexagram-dot', 318, '318',
                'hexagram-open-dot', 19, '19', 'star-triangle-up', 119,
                '119', 'star-triangle-up-open', 219, '219',
                'star-triangle-up-dot', 319, '319',
                'star-triangle-up-open-dot', 20, '20',
                'star-triangle-down', 120, '120',
                'star-triangle-down-open', 220, '220',
                'star-triangle-down-dot', 320, '320',
                'star-triangle-down-open-dot', 21, '21', 'star-square',
                121, '121', 'star-square-open', 221, '221',
                'star-square-dot', 321, '321', 'star-square-open-dot', 22,
                '22', 'star-diamond', 122, '122', 'star-diamond-open',
                222, '222', 'star-diamond-dot', 322, '322',
                'star-diamond-open-dot', 23, '23', 'diamond-tall', 123,
                '123', 'diamond-tall-open', 223, '223',
                'diamond-tall-dot', 323, '323', 'diamond-tall-open-dot',
                24, '24', 'diamond-wide', 124, '124', 'diamond-wide-open',
                224, '224', 'diamond-wide-dot', 324, '324',
                'diamond-wide-open-dot', 25, '25', 'hourglass', 125,
                '125', 'hourglass-open', 26, '26', 'bowtie', 126, '126',
                'bowtie-open', 27, '27', 'circle-cross', 127, '127',
                'circle-cross-open', 28, '28', 'circle-x', 128, '128',
                'circle-x-open', 29, '29', 'square-cross', 129, '129',
                'square-cross-open', 30, '30', 'square-x', 130, '130',
                'square-x-open', 31, '31', 'diamond-cross', 131, '131',
                'diamond-cross-open', 32, '32', 'diamond-x', 132, '132',
                'diamond-x-open', 33, '33', 'cross-thin', 133, '133',
                'cross-thin-open', 34, '34', 'x-thin', 134, '134',
                'x-thin-open', 35, '35', 'asterisk', 135, '135',
                'asterisk-open', 36, '36', 'hash', 136, '136',
                'hash-open', 236, '236', 'hash-dot', 336, '336',
                'hash-open-dot', 37, '37', 'y-up', 137, '137',
                'y-up-open', 38, '38', 'y-down', 138, '138',
                'y-down-open', 39, '39', 'y-left', 139, '139',
                'y-left-open', 40, '40', 'y-right', 140, '140',
                'y-right-open', 41, '41', 'line-ew', 141, '141',
                'line-ew-open', 42, '42', 'line-ns', 142, '142',
                'line-ns-open', 43, '43', 'line-ne', 143, '143',
                'line-ne-open', 44, '44', 'line-nw', 144, '144',
                'line-nw-open', 45, '45', 'arrow-up', 145, '145',
                'arrow-up-open', 46, '46', 'arrow-down', 146, '146',
                'arrow-down-open', 47, '47', 'arrow-left', 147, '147',
                'arrow-left-open', 48, '48', 'arrow-right', 148, '148',
                'arrow-right-open', 49, '49', 'arrow-bar-up', 149, '149',
                'arrow-bar-up-open', 50, '50', 'arrow-bar-down', 150,
                '150', 'arrow-bar-down-open', 51, '51', 'arrow-bar-left',
                151, '151', 'arrow-bar-left-open', 52, '52',
                'arrow-bar-right', 152, '152', 'arrow-bar-right-open', 53,
                '53', 'arrow', 153, '153', 'arrow-open', 54, '54',
                'arrow-wide', 154, '154', 'arrow-wide-open']
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["symbol"]

    @symbol.setter
    def symbol(self, val):
        self["symbol"] = val

    # symbolsrc
    # ---------
    @property
    def symbolsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `symbol`.

        The 'symbolsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["symbolsrc"]

    @symbolsrc.setter
    def symbolsrc(self, val):
        self["symbolsrc"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        angle
            Sets the marker angle in respect to `angleref`.
        angleref
            Sets the reference for marker angle. With "previous",
            angle 0 points along the line from the previous point
            to this one. With "up", angle 0 points toward the top
            of the screen.
        anglesrc
            Sets the source reference on Chart Studio Cloud for
            `angle`.
        autocolorscale
            Determines whether the colorscale is a default palette
            (`autocolorscale: true`) or the palette determined by
            `marker.colorscale`. Has an effect only if in
            `marker.color` is set to a numerical array. In case
            `colorscale` is unspecified or `autocolorscale` is
            true, the default palette will be chosen according to
            whether numbers in the `color` array are all positive,
            all negative or mixed.
        cauto
            Determines whether or not the color domain is computed
            with respect to the input data (here in `marker.color`)
            or the bounds set in `marker.cmin` and `marker.cmax`
            Has an effect only if in `marker.color` is set to a
            numerical array. Defaults to `false` when `marker.cmin`
            and `marker.cmax` are set by the user.
        cmax
            Sets the upper bound of the color domain. Has an effect
            only if in `marker.color` is set to a numerical array.
            Value should have the same units as in `marker.color`
            and if set, `marker.cmin` must be set as well.
        cmid
            Sets the mid-point of the color domain by scaling
            `marker.cmin` and/or `marker.cmax` to be equidistant to
            this point. Has an effect only if in `marker.color` is
            set to a numerical array. Value should have the same
            units as in `marker.color`. Has no effect when
            `marker.cauto` is `false`.
        cmin
            Sets the lower bound of the color domain. Has an effect
            only if in `marker.color` is set to a numerical array.
            Value should have the same units as in `marker.color`
            and if set, `marker.cmax` must be set as well.
        color
            Sets the marker color. It accepts either a specific
            color or an array of numbers that are mapped to the
            colorscale relative to the max and min values of the
            array or relative to `marker.cmin` and `marker.cmax` if
            set.
        coloraxis
            Sets a reference to a shared color axis. References to
            these shared color axes are "coloraxis", "coloraxis2",
            "coloraxis3", etc. Settings for these shared color axes
            are set in the layout, under `layout.coloraxis`,
            `layout.coloraxis2`, etc. Note that multiple color
            scales can be linked to the same color axis.
        colorbar
            :class:`plotly.graph_objects.scatter.marker.ColorBar`
            instance or dict with compatible properties
        colorscale
            Sets the colorscale. Has an effect only if in
            `marker.color` is set to a numerical array. The
            colorscale must be an array containing arrays mapping a
            normalized value to an rgb, rgba, hex, hsl, hsv, or
            named color string. At minimum, a mapping for the
            lowest (0) and highest (1) values are required. For
            example, `[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`.
            To control the bounds of the colorscale in color space,
            use `marker.cmin` and `marker.cmax`. Alternatively,
            `colorscale` may be a palette name string of the
            following list: Blackbody,Bluered,Blues,Cividis,Earth,E
            lectric,Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,Rd
            Bu,Reds,Viridis,YlGnBu,YlOrRd.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        gradient
            :class:`plotly.graph_objects.scatter.marker.Gradient`
            instance or dict with compatible properties
        line
            :class:`plotly.graph_objects.scatter.marker.Line`
            instance or dict with compatible properties
        maxdisplayed
            Sets a maximum number of points to be drawn on the
            graph. 0 corresponds to no limit.
        opacity
            Sets the marker opacity.
        opacitysrc
            Sets the source reference on Chart Studio Cloud for
            `opacity`.
        reversescale
            Reverses the color mapping if true. Has an effect only
            if in `marker.color` is set to a numerical array. If
            true, `marker.cmin` will correspond to the last color
            in the array and `marker.cmax` will correspond to the
            first color.
        showscale
            Determines whether or not a colorbar is displayed for
            this trace. Has an effect only if in `marker.color` is
            set to a numerical array.
        size
            Sets the marker size (in px).
        sizemin
            Has an effect only if `marker.size` is set to a
            numerical array. Sets the minimum size (in px) of the
            rendered marker points.
        sizemode
            Has an effect only if `marker.size` is set to a
            numerical array. Sets the rule for which the data in
            `size` is converted to pixels.
        sizeref
            Has an effect only if `marker.size` is set to a
            numerical array. Sets the scale factor used to
            determine the rendered size of marker points. Use with
            `sizemin` and `sizemode`.
        sizesrc
            Sets the source reference on Chart Studio Cloud for
            `size`.
        standoff
            Moves the marker away from the data point in the
            direction of `angle` (in px). This can be useful for
            example if you have another marker at this location and
            you want to point an arrowhead marker at it.
        standoffsrc
            Sets the source reference on Chart Studio Cloud for
            `standoff`.
        symbol
            Sets the marker symbol type. Adding 100 is equivalent
            to appending "-open" to a symbol name. Adding 200 is
            equivalent to appending "-dot" to a symbol name. Adding
            300 is equivalent to appending "-open-dot" or "dot-
            open" to a symbol name.
        symbolsrc
            Sets the source reference on Chart Studio Cloud for
            `symbol`.
        """

    def __init__(
        self,
        arg=None,
        angle=None,
        angleref=None,
        anglesrc=None,
        autocolorscale=None,
        cauto=None,
        cmax=None,
        cmid=None,
        cmin=None,
        color=None,
        coloraxis=None,
        colorbar=None,
        colorscale=None,
        colorsrc=None,
        gradient=None,
        line=None,
        maxdisplayed=None,
        opacity=None,
        opacitysrc=None,
        reversescale=None,
        showscale=None,
        size=None,
        sizemin=None,
        sizemode=None,
        sizeref=None,
        sizesrc=None,
        standoff=None,
        standoffsrc=None,
        symbol=None,
        symbolsrc=None,
        **kwargs,
    ):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.scatter.Marker`
        angle
            Sets the marker angle in respect to `angleref`.
        angleref
            Sets the reference for marker angle. With "previous",
            angle 0 points along the line from the previous point
            to this one. With "up", angle 0 points toward the top
            of the screen.
        anglesrc
            Sets the source reference on Chart Studio Cloud for
            `angle`.
        autocolorscale
            Determines whether the colorscale is a default palette
            (`autocolorscale: true`) or the palette determined by
            `marker.colorscale`. Has an effect only if in
            `marker.color` is set to a numerical array. In case
            `colorscale` is unspecified or `autocolorscale` is
            true, the default palette will be chosen according to
            whether numbers in the `color` array are all positive,
            all negative or mixed.
        cauto
            Determines whether or not the color domain is computed
            with respect to the input data (here in `marker.color`)
            or the bounds set in `marker.cmin` and `marker.cmax`
            Has an effect only if in `marker.color` is set to a
            numerical array. Defaults to `false` when `marker.cmin`
            and `marker.cmax` are set by the user.
        cmax
            Sets the upper bound of the color domain. Has an effect
            only if in `marker.color` is set to a numerical array.
            Value should have the same units as in `marker.color`
            and if set, `marker.cmin` must be set as well.
        cmid
            Sets the mid-point of the color domain by scaling
            `marker.cmin` and/or `marker.cmax` to be equidistant to
            this point. Has an effect only if in `marker.color` is
            set to a numerical array. Value should have the same
            units as in `marker.color`. Has no effect when
            `marker.cauto` is `false`.
        cmin
            Sets the lower bound of the color domain. Has an effect
            only if in `marker.color` is set to a numerical array.
            Value should have the same units as in `marker.color`
            and if set, `marker.cmax` must be set as well.
        color
            Sets the marker color. It accepts either a specific
            color or an array of numbers that are mapped to the
            colorscale relative to the max and min values of the
            array or relative to `marker.cmin` and `marker.cmax` if
            set.
        coloraxis
            Sets a reference to a shared color axis. References to
            these shared color axes are "coloraxis", "coloraxis2",
            "coloraxis3", etc. Settings for these shared color axes
            are set in the layout, under `layout.coloraxis`,
            `layout.coloraxis2`, etc. Note that multiple color
            scales can be linked to the same color axis.
        colorbar
            :class:`plotly.graph_objects.scatter.marker.ColorBar`
            instance or dict with compatible properties
        colorscale
            Sets the colorscale. Has an effect only if in
            `marker.color` is set to a numerical array. The
            colorscale must be an array containing arrays mapping a
            normalized value to an rgb, rgba, hex, hsl, hsv, or
            named color string. At minimum, a mapping for the
            lowest (0) and highest (1) values are required. For
            example, `[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`.
            To control the bounds of the colorscale in color space,
            use `marker.cmin` and `marker.cmax`. Alternatively,
            `colorscale` may be a palette name string of the
            following list: Blackbody,Bluered,Blues,Cividis,Earth,E
            lectric,Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,Rd
            Bu,Reds,Viridis,YlGnBu,YlOrRd.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        gradient
            :class:`plotly.graph_objects.scatter.marker.Gradient`
            instance or dict with compatible properties
        line
            :class:`plotly.graph_objects.scatter.marker.Line`
            instance or dict with compatible properties
        maxdisplayed
            Sets a maximum number of points to be drawn on the
            graph. 0 corresponds to no limit.
        opacity
            Sets the marker opacity.
        opacitysrc
            Sets the source reference on Chart Studio Cloud for
            `opacity`.
        reversescale
            Reverses the color mapping if true. Has an effect only
            if in `marker.color` is set to a numerical array. If
            true, `marker.cmin` will correspond to the last color
            in the array and `marker.cmax` will correspond to the
            first color.
        showscale
            Determines whether or not a colorbar is displayed for
            this trace. Has an effect only if in `marker.color` is
            set to a numerical array.
        size
            Sets the marker size (in px).
        sizemin
            Has an effect only if `marker.size` is set to a
            numerical array. Sets the minimum size (in px) of the
            rendered marker points.
        sizemode
            Has an effect only if `marker.size` is set to a
            numerical array. Sets the rule for which the data in
            `size` is converted to pixels.
        sizeref
            Has an effect only if `marker.size` is set to a
            numerical array. Sets the scale factor used to
            determine the rendered size of marker points. Use with
            `sizemin` and `sizemode`.
        sizesrc
            Sets the source reference on Chart Studio Cloud for
            `size`.
        standoff
            Moves the marker away from the data point in the
            direction of `angle` (in px). This can be useful for
            example if you have another marker at this location and
            you want to point an arrowhead marker at it.
        standoffsrc
            Sets the source reference on Chart Studio Cloud for
            `standoff`.
        symbol
            Sets the marker symbol type. Adding 100 is equivalent
            to appending "-open" to a symbol name. Adding 200 is
            equivalent to appending "-dot" to a symbol name. Adding
            300 is equivalent to appending "-open-dot" or "dot-
            open" to a symbol name.
        symbolsrc
            Sets the source reference on Chart Studio Cloud for
            `symbol`.

        Returns
        -------
        Marker
        """
        super(Marker, self).__init__("marker")

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
The first argument to the plotly.graph_objs.scatter.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scatter.Marker`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("angle", None)
        _v = angle if angle is not None else _v
        if _v is not None:
            self["angle"] = _v
        _v = arg.pop("angleref", None)
        _v = angleref if angleref is not None else _v
        if _v is not None:
            self["angleref"] = _v
        _v = arg.pop("anglesrc", None)
        _v = anglesrc if anglesrc is not None else _v
        if _v is not None:
            self["anglesrc"] = _v
        _v = arg.pop("autocolorscale", None)
        _v = autocolorscale if autocolorscale is not None else _v
        if _v is not None:
            self["autocolorscale"] = _v
        _v = arg.pop("cauto", None)
        _v = cauto if cauto is not None else _v
        if _v is not None:
            self["cauto"] = _v
        _v = arg.pop("cmax", None)
        _v = cmax if cmax is not None else _v
        if _v is not None:
            self["cmax"] = _v
        _v = arg.pop("cmid", None)
        _v = cmid if cmid is not None else _v
        if _v is not None:
            self["cmid"] = _v
        _v = arg.pop("cmin", None)
        _v = cmin if cmin is not None else _v
        if _v is not None:
            self["cmin"] = _v
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("coloraxis", None)
        _v = coloraxis if coloraxis is not None else _v
        if _v is not None:
            self["coloraxis"] = _v
        _v = arg.pop("colorbar", None)
        _v = colorbar if colorbar is not None else _v
        if _v is not None:
            self["colorbar"] = _v
        _v = arg.pop("colorscale", None)
        _v = colorscale if colorscale is not None else _v
        if _v is not None:
            self["colorscale"] = _v
        _v = arg.pop("colorsrc", None)
        _v = colorsrc if colorsrc is not None else _v
        if _v is not None:
            self["colorsrc"] = _v
        _v = arg.pop("gradient", None)
        _v = gradient if gradient is not None else _v
        if _v is not None:
            self["gradient"] = _v
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("maxdisplayed", None)
        _v = maxdisplayed if maxdisplayed is not None else _v
        if _v is not None:
            self["maxdisplayed"] = _v
        _v = arg.pop("opacity", None)
        _v = opacity if opacity is not None else _v
        if _v is not None:
            self["opacity"] = _v
        _v = arg.pop("opacitysrc", None)
        _v = opacitysrc if opacitysrc is not None else _v
        if _v is not None:
            self["opacitysrc"] = _v
        _v = arg.pop("reversescale", None)
        _v = reversescale if reversescale is not None else _v
        if _v is not None:
            self["reversescale"] = _v
        _v = arg.pop("showscale", None)
        _v = showscale if showscale is not None else _v
        if _v is not None:
            self["showscale"] = _v
        _v = arg.pop("size", None)
        _v = size if size is not None else _v
        if _v is not None:
            self["size"] = _v
        _v = arg.pop("sizemin", None)
        _v = sizemin if sizemin is not None else _v
        if _v is not None:
            self["sizemin"] = _v
        _v = arg.pop("sizemode", None)
        _v = sizemode if sizemode is not None else _v
        if _v is not None:
            self["sizemode"] = _v
        _v = arg.pop("sizeref", None)
        _v = sizeref if sizeref is not None else _v
        if _v is not None:
            self["sizeref"] = _v
        _v = arg.pop("sizesrc", None)
        _v = sizesrc if sizesrc is not None else _v
        if _v is not None:
            self["sizesrc"] = _v
        _v = arg.pop("standoff", None)
        _v = standoff if standoff is not None else _v
        if _v is not None:
            self["standoff"] = _v
        _v = arg.pop("standoffsrc", None)
        _v = standoffsrc if standoffsrc is not None else _v
        if _v is not None:
            self["standoffsrc"] = _v
        _v = arg.pop("symbol", None)
        _v = symbol if symbol is not None else _v
        if _v is not None:
            self["symbol"] = _v
        _v = arg.pop("symbolsrc", None)
        _v = symbolsrc if symbolsrc is not None else _v
        if _v is not None:
            self["symbolsrc"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
