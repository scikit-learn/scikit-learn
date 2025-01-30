from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "treemap"
    _path_str = "treemap.marker"
    _valid_props = {
        "autocolorscale",
        "cauto",
        "cmax",
        "cmid",
        "cmin",
        "coloraxis",
        "colorbar",
        "colors",
        "colorscale",
        "colorssrc",
        "cornerradius",
        "depthfade",
        "line",
        "pad",
        "pattern",
        "reversescale",
        "showscale",
    }

    # autocolorscale
    # --------------
    @property
    def autocolorscale(self):
        """
        Determines whether the colorscale is a default palette
        (`autocolorscale: true`) or the palette determined by
        `marker.colorscale`. Has an effect only if colors is set to a
        numerical array. In case `colorscale` is unspecified or
        `autocolorscale` is true, the default palette will be chosen
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
        respect to the input data (here colors) or the bounds set in
        `marker.cmin` and `marker.cmax` Has an effect only if colors is
        set to a numerical array. Defaults to `false` when
        `marker.cmin` and `marker.cmax` are set by the user.

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
        colors is set to a numerical array. Value should have the same
        units as colors and if set, `marker.cmin` must be set as well.

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
        effect only if colors is set to a numerical array. Value should
        have the same units as colors. Has no effect when
        `marker.cauto` is `false`.

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
        colors is set to a numerical array. Value should have the same
        units as colors and if set, `marker.cmax` must be set as well.

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
          - An instance of :class:`plotly.graph_objs.treemap.marker.ColorBar`
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
                    A tuple of :class:`plotly.graph_objects.treemap
                    .marker.colorbar.Tickformatstop` instances or
                    dicts with compatible properties
                tickformatstopdefaults
                    When used in a template (as layout.template.dat
                    a.treemap.marker.colorbar.tickformatstopdefault
                    s), sets the default property values to use for
                    elements of
                    treemap.marker.colorbar.tickformatstops
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
                    :class:`plotly.graph_objects.treemap.marker.col
                    orbar.Title` instance or dict with compatible
                    properties
                titlefont
                    Deprecated: Please use
                    treemap.marker.colorbar.title.font instead.
                    Sets this color bar's title font. Note that the
                    title's font used to be set by the now
                    deprecated `titlefont` attribute.
                titleside
                    Deprecated: Please use
                    treemap.marker.colorbar.title.side instead.
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
        plotly.graph_objs.treemap.marker.ColorBar
        """
        return self["colorbar"]

    @colorbar.setter
    def colorbar(self, val):
        self["colorbar"] = val

    # colors
    # ------
    @property
    def colors(self):
        """
        Sets the color of each sector of this trace. If not specified,
        the default trace color set is used to pick the sector colors.

        The 'colors' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["colors"]

    @colors.setter
    def colors(self, val):
        self["colors"] = val

    # colorscale
    # ----------
    @property
    def colorscale(self):
        """
        Sets the colorscale. Has an effect only if colors is set to a
        numerical array. The colorscale must be an array containing
        arrays mapping a normalized value to an rgb, rgba, hex, hsl,
        hsv, or named color string. At minimum, a mapping for the
        lowest (0) and highest (1) values are required. For example,
        `[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`. To control the
        bounds of the colorscale in color space, use `marker.cmin` and
        `marker.cmax`. Alternatively, `colorscale` may be a palette
        name string of the following list: Blackbody,Bluered,Blues,Civi
        dis,Earth,Electric,Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow
        ,RdBu,Reds,Viridis,YlGnBu,YlOrRd.

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

    # colorssrc
    # ---------
    @property
    def colorssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `colors`.

        The 'colorssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["colorssrc"]

    @colorssrc.setter
    def colorssrc(self, val):
        self["colorssrc"] = val

    # cornerradius
    # ------------
    @property
    def cornerradius(self):
        """
        Sets the maximum rounding of corners (in px).

        The 'cornerradius' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["cornerradius"]

    @cornerradius.setter
    def cornerradius(self, val):
        self["cornerradius"] = val

    # depthfade
    # ---------
    @property
    def depthfade(self):
        """
        Determines if the sector colors are faded towards the
        background from the leaves up to the headers. This option is
        unavailable when a `colorscale` is present, defaults to false
        when `marker.colors` is set, but otherwise defaults to true.
        When set to "reversed", the fading direction is inverted, that
        is the top elements within hierarchy are drawn with fully
        saturated colors while the leaves are faded towards the
        background color.

        The 'depthfade' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, False, 'reversed']

        Returns
        -------
        Any
        """
        return self["depthfade"]

    @depthfade.setter
    def depthfade(self, val):
        self["depthfade"] = val

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.treemap.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the color of the line enclosing each
                    sector. Defaults to the `paper_bgcolor` value.
                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.
                width
                    Sets the width (in px) of the line enclosing
                    each sector.
                widthsrc
                    Sets the source reference on Chart Studio Cloud
                    for `width`.

        Returns
        -------
        plotly.graph_objs.treemap.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # pad
    # ---
    @property
    def pad(self):
        """
        The 'pad' property is an instance of Pad
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.treemap.marker.Pad`
          - A dict of string/value properties that will be passed
            to the Pad constructor

            Supported dict properties:

                b
                    Sets the padding form the bottom (in px).
                l
                    Sets the padding form the left (in px).
                r
                    Sets the padding form the right (in px).
                t
                    Sets the padding form the top (in px).

        Returns
        -------
        plotly.graph_objs.treemap.marker.Pad
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    # pattern
    # -------
    @property
    def pattern(self):
        """
        Sets the pattern within the marker.

        The 'pattern' property is an instance of Pattern
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.treemap.marker.Pattern`
          - A dict of string/value properties that will be passed
            to the Pattern constructor

            Supported dict properties:

                bgcolor
                    When there is no colorscale sets the color of
                    background pattern fill. Defaults to a
                    `marker.color` background when `fillmode` is
                    "overlay". Otherwise, defaults to a transparent
                    background.
                bgcolorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `bgcolor`.
                fgcolor
                    When there is no colorscale sets the color of
                    foreground pattern fill. Defaults to a
                    `marker.color` background when `fillmode` is
                    "replace". Otherwise, defaults to dark grey or
                    white to increase contrast with the `bgcolor`.
                fgcolorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `fgcolor`.
                fgopacity
                    Sets the opacity of the foreground pattern
                    fill. Defaults to a 0.5 when `fillmode` is
                    "overlay". Otherwise, defaults to 1.
                fillmode
                    Determines whether `marker.color` should be
                    used as a default to `bgcolor` or a `fgcolor`.
                shape
                    Sets the shape of the pattern fill. By default,
                    no pattern is used for filling the area.
                shapesrc
                    Sets the source reference on Chart Studio Cloud
                    for `shape`.
                size
                    Sets the size of unit squares of the pattern
                    fill in pixels, which corresponds to the
                    interval of repetition of the pattern.
                sizesrc
                    Sets the source reference on Chart Studio Cloud
                    for `size`.
                solidity
                    Sets the solidity of the pattern fill. Solidity
                    is roughly the fraction of the area filled by
                    the pattern. Solidity of 0 shows only the
                    background color without pattern and solidty of
                    1 shows only the foreground color without
                    pattern.
                soliditysrc
                    Sets the source reference on Chart Studio Cloud
                    for `solidity`.

        Returns
        -------
        plotly.graph_objs.treemap.marker.Pattern
        """
        return self["pattern"]

    @pattern.setter
    def pattern(self, val):
        self["pattern"] = val

    # reversescale
    # ------------
    @property
    def reversescale(self):
        """
        Reverses the color mapping if true. Has an effect only if
        colors is set to a numerical array. If true, `marker.cmin` will
        correspond to the last color in the array and `marker.cmax`
        will correspond to the first color.

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
        trace. Has an effect only if colors is set to a numerical
        array.

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

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        autocolorscale
            Determines whether the colorscale is a default palette
            (`autocolorscale: true`) or the palette determined by
            `marker.colorscale`. Has an effect only if colors is
            set to a numerical array. In case `colorscale` is
            unspecified or `autocolorscale` is true, the default
            palette will be chosen according to whether numbers in
            the `color` array are all positive, all negative or
            mixed.
        cauto
            Determines whether or not the color domain is computed
            with respect to the input data (here colors) or the
            bounds set in `marker.cmin` and `marker.cmax` Has an
            effect only if colors is set to a numerical array.
            Defaults to `false` when `marker.cmin` and
            `marker.cmax` are set by the user.
        cmax
            Sets the upper bound of the color domain. Has an effect
            only if colors is set to a numerical array. Value
            should have the same units as colors and if set,
            `marker.cmin` must be set as well.
        cmid
            Sets the mid-point of the color domain by scaling
            `marker.cmin` and/or `marker.cmax` to be equidistant to
            this point. Has an effect only if colors is set to a
            numerical array. Value should have the same units as
            colors. Has no effect when `marker.cauto` is `false`.
        cmin
            Sets the lower bound of the color domain. Has an effect
            only if colors is set to a numerical array. Value
            should have the same units as colors and if set,
            `marker.cmax` must be set as well.
        coloraxis
            Sets a reference to a shared color axis. References to
            these shared color axes are "coloraxis", "coloraxis2",
            "coloraxis3", etc. Settings for these shared color axes
            are set in the layout, under `layout.coloraxis`,
            `layout.coloraxis2`, etc. Note that multiple color
            scales can be linked to the same color axis.
        colorbar
            :class:`plotly.graph_objects.treemap.marker.ColorBar`
            instance or dict with compatible properties
        colors
            Sets the color of each sector of this trace. If not
            specified, the default trace color set is used to pick
            the sector colors.
        colorscale
            Sets the colorscale. Has an effect only if colors is
            set to a numerical array. The colorscale must be an
            array containing arrays mapping a normalized value to
            an rgb, rgba, hex, hsl, hsv, or named color string. At
            minimum, a mapping for the lowest (0) and highest (1)
            values are required. For example, `[[0,
            'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`. To control the
            bounds of the colorscale in color space, use
            `marker.cmin` and `marker.cmax`. Alternatively,
            `colorscale` may be a palette name string of the
            following list: Blackbody,Bluered,Blues,Cividis,Earth,E
            lectric,Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,Rd
            Bu,Reds,Viridis,YlGnBu,YlOrRd.
        colorssrc
            Sets the source reference on Chart Studio Cloud for
            `colors`.
        cornerradius
            Sets the maximum rounding of corners (in px).
        depthfade
            Determines if the sector colors are faded towards the
            background from the leaves up to the headers. This
            option is unavailable when a `colorscale` is present,
            defaults to false when `marker.colors` is set, but
            otherwise defaults to true. When set to "reversed", the
            fading direction is inverted, that is the top elements
            within hierarchy are drawn with fully saturated colors
            while the leaves are faded towards the background
            color.
        line
            :class:`plotly.graph_objects.treemap.marker.Line`
            instance or dict with compatible properties
        pad
            :class:`plotly.graph_objects.treemap.marker.Pad`
            instance or dict with compatible properties
        pattern
            Sets the pattern within the marker.
        reversescale
            Reverses the color mapping if true. Has an effect only
            if colors is set to a numerical array. If true,
            `marker.cmin` will correspond to the last color in the
            array and `marker.cmax` will correspond to the first
            color.
        showscale
            Determines whether or not a colorbar is displayed for
            this trace. Has an effect only if colors is set to a
            numerical array.
        """

    def __init__(
        self,
        arg=None,
        autocolorscale=None,
        cauto=None,
        cmax=None,
        cmid=None,
        cmin=None,
        coloraxis=None,
        colorbar=None,
        colors=None,
        colorscale=None,
        colorssrc=None,
        cornerradius=None,
        depthfade=None,
        line=None,
        pad=None,
        pattern=None,
        reversescale=None,
        showscale=None,
        **kwargs,
    ):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.treemap.Marker`
        autocolorscale
            Determines whether the colorscale is a default palette
            (`autocolorscale: true`) or the palette determined by
            `marker.colorscale`. Has an effect only if colors is
            set to a numerical array. In case `colorscale` is
            unspecified or `autocolorscale` is true, the default
            palette will be chosen according to whether numbers in
            the `color` array are all positive, all negative or
            mixed.
        cauto
            Determines whether or not the color domain is computed
            with respect to the input data (here colors) or the
            bounds set in `marker.cmin` and `marker.cmax` Has an
            effect only if colors is set to a numerical array.
            Defaults to `false` when `marker.cmin` and
            `marker.cmax` are set by the user.
        cmax
            Sets the upper bound of the color domain. Has an effect
            only if colors is set to a numerical array. Value
            should have the same units as colors and if set,
            `marker.cmin` must be set as well.
        cmid
            Sets the mid-point of the color domain by scaling
            `marker.cmin` and/or `marker.cmax` to be equidistant to
            this point. Has an effect only if colors is set to a
            numerical array. Value should have the same units as
            colors. Has no effect when `marker.cauto` is `false`.
        cmin
            Sets the lower bound of the color domain. Has an effect
            only if colors is set to a numerical array. Value
            should have the same units as colors and if set,
            `marker.cmax` must be set as well.
        coloraxis
            Sets a reference to a shared color axis. References to
            these shared color axes are "coloraxis", "coloraxis2",
            "coloraxis3", etc. Settings for these shared color axes
            are set in the layout, under `layout.coloraxis`,
            `layout.coloraxis2`, etc. Note that multiple color
            scales can be linked to the same color axis.
        colorbar
            :class:`plotly.graph_objects.treemap.marker.ColorBar`
            instance or dict with compatible properties
        colors
            Sets the color of each sector of this trace. If not
            specified, the default trace color set is used to pick
            the sector colors.
        colorscale
            Sets the colorscale. Has an effect only if colors is
            set to a numerical array. The colorscale must be an
            array containing arrays mapping a normalized value to
            an rgb, rgba, hex, hsl, hsv, or named color string. At
            minimum, a mapping for the lowest (0) and highest (1)
            values are required. For example, `[[0,
            'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`. To control the
            bounds of the colorscale in color space, use
            `marker.cmin` and `marker.cmax`. Alternatively,
            `colorscale` may be a palette name string of the
            following list: Blackbody,Bluered,Blues,Cividis,Earth,E
            lectric,Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,Rd
            Bu,Reds,Viridis,YlGnBu,YlOrRd.
        colorssrc
            Sets the source reference on Chart Studio Cloud for
            `colors`.
        cornerradius
            Sets the maximum rounding of corners (in px).
        depthfade
            Determines if the sector colors are faded towards the
            background from the leaves up to the headers. This
            option is unavailable when a `colorscale` is present,
            defaults to false when `marker.colors` is set, but
            otherwise defaults to true. When set to "reversed", the
            fading direction is inverted, that is the top elements
            within hierarchy are drawn with fully saturated colors
            while the leaves are faded towards the background
            color.
        line
            :class:`plotly.graph_objects.treemap.marker.Line`
            instance or dict with compatible properties
        pad
            :class:`plotly.graph_objects.treemap.marker.Pad`
            instance or dict with compatible properties
        pattern
            Sets the pattern within the marker.
        reversescale
            Reverses the color mapping if true. Has an effect only
            if colors is set to a numerical array. If true,
            `marker.cmin` will correspond to the last color in the
            array and `marker.cmax` will correspond to the first
            color.
        showscale
            Determines whether or not a colorbar is displayed for
            this trace. Has an effect only if colors is set to a
            numerical array.

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
The first argument to the plotly.graph_objs.treemap.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.treemap.Marker`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
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
        _v = arg.pop("coloraxis", None)
        _v = coloraxis if coloraxis is not None else _v
        if _v is not None:
            self["coloraxis"] = _v
        _v = arg.pop("colorbar", None)
        _v = colorbar if colorbar is not None else _v
        if _v is not None:
            self["colorbar"] = _v
        _v = arg.pop("colors", None)
        _v = colors if colors is not None else _v
        if _v is not None:
            self["colors"] = _v
        _v = arg.pop("colorscale", None)
        _v = colorscale if colorscale is not None else _v
        if _v is not None:
            self["colorscale"] = _v
        _v = arg.pop("colorssrc", None)
        _v = colorssrc if colorssrc is not None else _v
        if _v is not None:
            self["colorssrc"] = _v
        _v = arg.pop("cornerradius", None)
        _v = cornerradius if cornerradius is not None else _v
        if _v is not None:
            self["cornerradius"] = _v
        _v = arg.pop("depthfade", None)
        _v = depthfade if depthfade is not None else _v
        if _v is not None:
            self["depthfade"] = _v
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("pad", None)
        _v = pad if pad is not None else _v
        if _v is not None:
            self["pad"] = _v
        _v = arg.pop("pattern", None)
        _v = pattern if pattern is not None else _v
        if _v is not None:
            self["pattern"] = _v
        _v = arg.pop("reversescale", None)
        _v = reversescale if reversescale is not None else _v
        if _v is not None:
            self["reversescale"] = _v
        _v = arg.pop("showscale", None)
        _v = showscale if showscale is not None else _v
        if _v is not None:
            self["showscale"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
