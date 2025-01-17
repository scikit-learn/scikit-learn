from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class XAxis(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.xaxis"
    _valid_props = {
        "anchor",
        "automargin",
        "autorange",
        "autorangeoptions",
        "autotickangles",
        "autotypenumbers",
        "calendar",
        "categoryarray",
        "categoryarraysrc",
        "categoryorder",
        "color",
        "constrain",
        "constraintoward",
        "dividercolor",
        "dividerwidth",
        "domain",
        "dtick",
        "exponentformat",
        "fixedrange",
        "gridcolor",
        "griddash",
        "gridwidth",
        "hoverformat",
        "insiderange",
        "labelalias",
        "layer",
        "linecolor",
        "linewidth",
        "matches",
        "maxallowed",
        "minallowed",
        "minexponent",
        "minor",
        "mirror",
        "nticks",
        "overlaying",
        "position",
        "range",
        "rangebreakdefaults",
        "rangebreaks",
        "rangemode",
        "rangeselector",
        "rangeslider",
        "scaleanchor",
        "scaleratio",
        "separatethousands",
        "showdividers",
        "showexponent",
        "showgrid",
        "showline",
        "showspikes",
        "showticklabels",
        "showtickprefix",
        "showticksuffix",
        "side",
        "spikecolor",
        "spikedash",
        "spikemode",
        "spikesnap",
        "spikethickness",
        "tick0",
        "tickangle",
        "tickcolor",
        "tickfont",
        "tickformat",
        "tickformatstopdefaults",
        "tickformatstops",
        "ticklabelindex",
        "ticklabelindexsrc",
        "ticklabelmode",
        "ticklabeloverflow",
        "ticklabelposition",
        "ticklabelshift",
        "ticklabelstandoff",
        "ticklabelstep",
        "ticklen",
        "tickmode",
        "tickprefix",
        "ticks",
        "tickson",
        "ticksuffix",
        "ticktext",
        "ticktextsrc",
        "tickvals",
        "tickvalssrc",
        "tickwidth",
        "title",
        "titlefont",
        "type",
        "uirevision",
        "visible",
        "zeroline",
        "zerolinecolor",
        "zerolinewidth",
    }

    # anchor
    # ------
    @property
    def anchor(self):
        """
        If set to an opposite-letter axis id (e.g. `x2`, `y`), this
        axis is bound to the corresponding opposite-letter axis. If set
        to "free", this axis' position is determined by `position`.

        The 'anchor' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['free']
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$',
                '^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["anchor"]

    @anchor.setter
    def anchor(self, val):
        self["anchor"] = val

    # automargin
    # ----------
    @property
    def automargin(self):
        """
        Determines whether long tick labels automatically grow the
        figure margins.

        The 'automargin' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['height', 'width', 'left', 'right', 'top', 'bottom'] joined with '+' characters
            (e.g. 'height+width')
            OR exactly one of [True, False] (e.g. 'False')

        Returns
        -------
        Any
        """
        return self["automargin"]

    @automargin.setter
    def automargin(self, val):
        self["automargin"] = val

    # autorange
    # ---------
    @property
    def autorange(self):
        """
        Determines whether or not the range of this axis is computed in
        relation to the input data. See `rangemode` for more info. If
        `range` is provided and it has a value for both the lower and
        upper bound, `autorange` is set to False. Using "min" applies
        autorange only to set the minimum. Using "max" applies
        autorange only to set the maximum. Using *min reversed* applies
        autorange only to set the minimum on a reversed axis. Using
        *max reversed* applies autorange only to set the maximum on a
        reversed axis. Using "reversed" applies autorange on both ends
        and reverses the axis direction.

        The 'autorange' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, False, 'reversed', 'min reversed', 'max reversed',
                'min', 'max']

        Returns
        -------
        Any
        """
        return self["autorange"]

    @autorange.setter
    def autorange(self, val):
        self["autorange"] = val

    # autorangeoptions
    # ----------------
    @property
    def autorangeoptions(self):
        """
        The 'autorangeoptions' property is an instance of Autorangeoptions
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Autorangeoptions`
          - A dict of string/value properties that will be passed
            to the Autorangeoptions constructor

            Supported dict properties:

                clipmax
                    Clip autorange maximum if it goes beyond this
                    value. Has no effect when
                    `autorangeoptions.maxallowed` is provided.
                clipmin
                    Clip autorange minimum if it goes beyond this
                    value. Has no effect when
                    `autorangeoptions.minallowed` is provided.
                include
                    Ensure this value is included in autorange.
                includesrc
                    Sets the source reference on Chart Studio Cloud
                    for `include`.
                maxallowed
                    Use this value exactly as autorange maximum.
                minallowed
                    Use this value exactly as autorange minimum.

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Autorangeoptions
        """
        return self["autorangeoptions"]

    @autorangeoptions.setter
    def autorangeoptions(self, val):
        self["autorangeoptions"] = val

    # autotickangles
    # --------------
    @property
    def autotickangles(self):
        """
        When `tickangle` is set to "auto", it will be set to the first
        angle in this array that is large enough to prevent label
        overlap.

        The 'autotickangles' property is an info array that may be specified as:
        * a list of elements where:
          The 'autotickangles[i]' property is a angle (in degrees) that may be
        specified as a number between -180 and 180.
        Numeric values outside this range are converted to the equivalent value
        (e.g. 270 is converted to -90).

        Returns
        -------
        list
        """
        return self["autotickangles"]

    @autotickangles.setter
    def autotickangles(self, val):
        self["autotickangles"] = val

    # autotypenumbers
    # ---------------
    @property
    def autotypenumbers(self):
        """
        Using "strict" a numeric string in trace data is not converted
        to a number. Using *convert types* a numeric string in trace
        data may be treated as a number during automatic axis `type`
        detection. Defaults to layout.autotypenumbers.

        The 'autotypenumbers' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['convert types', 'strict']

        Returns
        -------
        Any
        """
        return self["autotypenumbers"]

    @autotypenumbers.setter
    def autotypenumbers(self, val):
        self["autotypenumbers"] = val

    # calendar
    # --------
    @property
    def calendar(self):
        """
        Sets the calendar system to use for `range` and `tick0` if this
        is a date axis. This does not set the calendar for interpreting
        data on this axis, that's specified in the trace or via the
        global `layout.calendar`

        The 'calendar' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['chinese', 'coptic', 'discworld', 'ethiopian',
                'gregorian', 'hebrew', 'islamic', 'jalali', 'julian',
                'mayan', 'nanakshahi', 'nepali', 'persian', 'taiwan',
                'thai', 'ummalqura']

        Returns
        -------
        Any
        """
        return self["calendar"]

    @calendar.setter
    def calendar(self, val):
        self["calendar"] = val

    # categoryarray
    # -------------
    @property
    def categoryarray(self):
        """
        Sets the order in which categories on this axis appear. Only
        has an effect if `categoryorder` is set to "array". Used with
        `categoryorder`.

        The 'categoryarray' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["categoryarray"]

    @categoryarray.setter
    def categoryarray(self, val):
        self["categoryarray"] = val

    # categoryarraysrc
    # ----------------
    @property
    def categoryarraysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `categoryarray`.

        The 'categoryarraysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["categoryarraysrc"]

    @categoryarraysrc.setter
    def categoryarraysrc(self, val):
        self["categoryarraysrc"] = val

    # categoryorder
    # -------------
    @property
    def categoryorder(self):
        """
        Specifies the ordering logic for the case of categorical
        variables. By default, plotly uses "trace", which specifies the
        order that is present in the data supplied. Set `categoryorder`
        to *category ascending* or *category descending* if order
        should be determined by the alphanumerical order of the
        category names. Set `categoryorder` to "array" to derive the
        ordering from the attribute `categoryarray`. If a category is
        not found in the `categoryarray` array, the sorting behavior
        for that attribute will be identical to the "trace" mode. The
        unspecified categories will follow the categories in
        `categoryarray`. Set `categoryorder` to *total ascending* or
        *total descending* if order should be determined by the
        numerical order of the values. Similarly, the order can be
        determined by the min, max, sum, mean, geometric mean or median
        of all the values.

        The 'categoryorder' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['trace', 'category ascending', 'category descending',
                'array', 'total ascending', 'total descending', 'min
                ascending', 'min descending', 'max ascending', 'max
                descending', 'sum ascending', 'sum descending', 'mean
                ascending', 'mean descending', 'geometric mean ascending',
                'geometric mean descending', 'median ascending', 'median
                descending']

        Returns
        -------
        Any
        """
        return self["categoryorder"]

    @categoryorder.setter
    def categoryorder(self, val):
        self["categoryorder"] = val

    # color
    # -----
    @property
    def color(self):
        """
        Sets default for all colors associated with this axis all at
        once: line, font, tick, and grid colors. Grid color is
        lightened by blending this with the plot background Individual
        pieces can override this.

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

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    # constrain
    # ---------
    @property
    def constrain(self):
        """
        If this axis needs to be compressed (either due to its own
        `scaleanchor` and `scaleratio` or those of the other axis),
        determines how that happens: by increasing the "range", or by
        decreasing the "domain". Default is "domain" for axes
        containing image traces, "range" otherwise.

        The 'constrain' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['range', 'domain']

        Returns
        -------
        Any
        """
        return self["constrain"]

    @constrain.setter
    def constrain(self, val):
        self["constrain"] = val

    # constraintoward
    # ---------------
    @property
    def constraintoward(self):
        """
        If this axis needs to be compressed (either due to its own
        `scaleanchor` and `scaleratio` or those of the other axis),
        determines which direction we push the originally specified
        plot area. Options are "left", "center" (default), and "right"
        for x axes, and "top", "middle" (default), and "bottom" for y
        axes.

        The 'constraintoward' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['left', 'center', 'right', 'top', 'middle', 'bottom']

        Returns
        -------
        Any
        """
        return self["constraintoward"]

    @constraintoward.setter
    def constraintoward(self, val):
        self["constraintoward"] = val

    # dividercolor
    # ------------
    @property
    def dividercolor(self):
        """
        Sets the color of the dividers Only has an effect on
        "multicategory" axes.

        The 'dividercolor' property is a color and may be specified as:
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
        return self["dividercolor"]

    @dividercolor.setter
    def dividercolor(self, val):
        self["dividercolor"] = val

    # dividerwidth
    # ------------
    @property
    def dividerwidth(self):
        """
        Sets the width (in px) of the dividers Only has an effect on
        "multicategory" axes.

        The 'dividerwidth' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["dividerwidth"]

    @dividerwidth.setter
    def dividerwidth(self, val):
        self["dividerwidth"] = val

    # domain
    # ------
    @property
    def domain(self):
        """
            Sets the domain of this axis (in plot fraction).

            The 'domain' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'domain[0]' property is a number and may be specified as:
              - An int or float in the interval [0, 1]
        (1) The 'domain[1]' property is a number and may be specified as:
              - An int or float in the interval [0, 1]

            Returns
            -------
            list
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    # dtick
    # -----
    @property
    def dtick(self):
        """
        Sets the step in-between ticks on this axis. Use with `tick0`.
        Must be a positive number, or special strings available to
        "log" and "date" axes. If the axis `type` is "log", then ticks
        are set every 10^(n*dtick) where n is the tick number. For
        example, to set a tick mark at 1, 10, 100, 1000, ... set dtick
        to 1. To set tick marks at 1, 100, 10000, ... set dtick to 2.
        To set tick marks at 1, 5, 25, 125, 625, 3125, ... set dtick to
        log_10(5), or 0.69897000433. "log" has several special values;
        "L<f>", where `f` is a positive number, gives ticks linearly
        spaced in value (but not position). For example `tick0` = 0.1,
        `dtick` = "L0.5" will put ticks at 0.1, 0.6, 1.1, 1.6 etc. To
        show powers of 10 plus small digits between, use "D1" (all
        digits) or "D2" (only 2 and 5). `tick0` is ignored for "D1" and
        "D2". If the axis `type` is "date", then you must convert the
        time to milliseconds. For example, to set the interval between
        ticks to one day, set `dtick` to 86400000.0. "date" also has
        special values "M<n>" gives ticks spaced by a number of months.
        `n` must be a positive integer. To set ticks on the 15th of
        every third month, set `tick0` to "2000-01-15" and `dtick` to
        "M3". To set ticks every 4 years, set `dtick` to "M48"

        The 'dtick' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["dtick"]

    @dtick.setter
    def dtick(self, val):
        self["dtick"] = val

    # exponentformat
    # --------------
    @property
    def exponentformat(self):
        """
        Determines a formatting rule for the tick exponents. For
        example, consider the number 1,000,000,000. If "none", it
        appears as 1,000,000,000. If "e", 1e+9. If "E", 1E+9. If
        "power", 1x10^9 (with 9 in a super script). If "SI", 1G. If
        "B", 1B.

        The 'exponentformat' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['none', 'e', 'E', 'power', 'SI', 'B']

        Returns
        -------
        Any
        """
        return self["exponentformat"]

    @exponentformat.setter
    def exponentformat(self, val):
        self["exponentformat"] = val

    # fixedrange
    # ----------
    @property
    def fixedrange(self):
        """
        Determines whether or not this axis is zoom-able. If true, then
        zoom is disabled.

        The 'fixedrange' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["fixedrange"]

    @fixedrange.setter
    def fixedrange(self, val):
        self["fixedrange"] = val

    # gridcolor
    # ---------
    @property
    def gridcolor(self):
        """
        Sets the color of the grid lines.

        The 'gridcolor' property is a color and may be specified as:
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
        return self["gridcolor"]

    @gridcolor.setter
    def gridcolor(self, val):
        self["gridcolor"] = val

    # griddash
    # --------
    @property
    def griddash(self):
        """
        Sets the dash style of lines. Set to a dash type string
        ("solid", "dot", "dash", "longdash", "dashdot", or
        "longdashdot") or a dash length list in px (eg
        "5px,10px,2px,2px").

        The 'griddash' property is an enumeration that may be specified as:
          - One of the following dash styles:
                ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
          - A string containing a dash length list in pixels or percentages
                (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)

        Returns
        -------
        str
        """
        return self["griddash"]

    @griddash.setter
    def griddash(self, val):
        self["griddash"] = val

    # gridwidth
    # ---------
    @property
    def gridwidth(self):
        """
        Sets the width (in px) of the grid lines.

        The 'gridwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["gridwidth"]

    @gridwidth.setter
    def gridwidth(self, val):
        self["gridwidth"] = val

    # hoverformat
    # -----------
    @property
    def hoverformat(self):
        """
        Sets the hover text formatting rule using d3 formatting mini-
        languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format. And for
        dates see: https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format. We add two items to d3's date
        formatter: "%h" for half of the year as a decimal number as
        well as "%{n}f" for fractional seconds with n digits. For
        example, *2016-10-13 09:15:23.456* with tickformat
        "%H~%M~%S.%2f" would display "09~15~23.46"

        The 'hoverformat' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["hoverformat"]

    @hoverformat.setter
    def hoverformat(self, val):
        self["hoverformat"] = val

    # insiderange
    # -----------
    @property
    def insiderange(self):
        """
            Could be used to set the desired inside range of this axis
            (excluding the labels) when `ticklabelposition` of the anchored
            axis has "inside". Not implemented for axes with `type` "log".
            This would be ignored when `range` is provided.

            The 'insiderange' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'insiderange[0]' property accepts values of any type
        (1) The 'insiderange[1]' property accepts values of any type

            Returns
            -------
            list
        """
        return self["insiderange"]

    @insiderange.setter
    def insiderange(self, val):
        self["insiderange"] = val

    # labelalias
    # ----------
    @property
    def labelalias(self):
        """
        Replacement text for specific tick or hover labels. For example
        using {US: 'USA', CA: 'Canada'} changes US to USA and CA to
        Canada. The labels we would have shown must match the keys
        exactly, after adding any tickprefix or ticksuffix. For
        negative numbers the minus sign symbol used (U+2212) is wider
        than the regular ascii dash. That means you need to use −1
        instead of -1. labelalias can be used with any axis type, and
        both keys (if needed) and values (if desired) can include html-
        like tags or MathJax.

        The 'labelalias' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["labelalias"]

    @labelalias.setter
    def labelalias(self, val):
        self["labelalias"] = val

    # layer
    # -----
    @property
    def layer(self):
        """
        Sets the layer on which this axis is displayed. If *above
        traces*, this axis is displayed above all the subplot's traces
        If *below traces*, this axis is displayed below all the
        subplot's traces, but above the grid lines. Useful when used
        together with scatter-like traces with `cliponaxis` set to
        False to show markers and/or text nodes above this axis.

        The 'layer' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['above traces', 'below traces']

        Returns
        -------
        Any
        """
        return self["layer"]

    @layer.setter
    def layer(self, val):
        self["layer"] = val

    # linecolor
    # ---------
    @property
    def linecolor(self):
        """
        Sets the axis line color.

        The 'linecolor' property is a color and may be specified as:
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
        return self["linecolor"]

    @linecolor.setter
    def linecolor(self, val):
        self["linecolor"] = val

    # linewidth
    # ---------
    @property
    def linewidth(self):
        """
        Sets the width (in px) of the axis line.

        The 'linewidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["linewidth"]

    @linewidth.setter
    def linewidth(self, val):
        self["linewidth"] = val

    # matches
    # -------
    @property
    def matches(self):
        """
        If set to another axis id (e.g. `x2`, `y`), the range of this
        axis will match the range of the corresponding axis in data-
        coordinates space. Moreover, matching axes share auto-range
        values, category lists and histogram auto-bins. Note that
        setting axes simultaneously in both a `scaleanchor` and a
        `matches` constraint is currently forbidden. Moreover, note
        that matching axes must have the same `type`.

        The 'matches' property is an enumeration that may be specified as:
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$',
                '^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["matches"]

    @matches.setter
    def matches(self, val):
        self["matches"] = val

    # maxallowed
    # ----------
    @property
    def maxallowed(self):
        """
        Determines the maximum range of this axis.

        The 'maxallowed' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["maxallowed"]

    @maxallowed.setter
    def maxallowed(self, val):
        self["maxallowed"] = val

    # minallowed
    # ----------
    @property
    def minallowed(self):
        """
        Determines the minimum range of this axis.

        The 'minallowed' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["minallowed"]

    @minallowed.setter
    def minallowed(self, val):
        self["minallowed"] = val

    # minexponent
    # -----------
    @property
    def minexponent(self):
        """
        Hide SI prefix for 10^n if |n| is below this number. This only
        has an effect when `tickformat` is "SI" or "B".

        The 'minexponent' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["minexponent"]

    @minexponent.setter
    def minexponent(self, val):
        self["minexponent"] = val

    # minor
    # -----
    @property
    def minor(self):
        """
        The 'minor' property is an instance of Minor
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Minor`
          - A dict of string/value properties that will be passed
            to the Minor constructor

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
                gridcolor
                    Sets the color of the grid lines.
                griddash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                gridwidth
                    Sets the width (in px) of the grid lines.
                nticks
                    Specifies the maximum number of ticks for the
                    particular axis. The actual number of ticks
                    will be chosen automatically to be less than or
                    equal to `nticks`. Has an effect only if
                    `tickmode` is set to "auto".
                showgrid
                    Determines whether or not grid lines are drawn.
                    If True, the grid lines are drawn at every tick
                    mark.
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
                tickcolor
                    Sets the tick color.
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
                ticks
                    Determines whether ticks are drawn or not. If
                    "", this axis' ticks are not drawn. If
                    "outside" ("inside"), this axis' are drawn
                    outside (inside) the axis lines.
                tickvals
                    Sets the values at which ticks on this axis
                    appear. Only has an effect if `tickmode` is set
                    to "array". Used with `ticktext`.
                tickvalssrc
                    Sets the source reference on Chart Studio Cloud
                    for `tickvals`.
                tickwidth
                    Sets the tick width (in px).

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Minor
        """
        return self["minor"]

    @minor.setter
    def minor(self, val):
        self["minor"] = val

    # mirror
    # ------
    @property
    def mirror(self):
        """
        Determines if the axis lines or/and ticks are mirrored to the
        opposite side of the plotting area. If True, the axis lines are
        mirrored. If "ticks", the axis lines and ticks are mirrored. If
        False, mirroring is disable. If "all", axis lines are mirrored
        on all shared-axes subplots. If "allticks", axis lines and
        ticks are mirrored on all shared-axes subplots.

        The 'mirror' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, 'ticks', False, 'all', 'allticks']

        Returns
        -------
        Any
        """
        return self["mirror"]

    @mirror.setter
    def mirror(self, val):
        self["mirror"] = val

    # nticks
    # ------
    @property
    def nticks(self):
        """
        Specifies the maximum number of ticks for the particular axis.
        The actual number of ticks will be chosen automatically to be
        less than or equal to `nticks`. Has an effect only if
        `tickmode` is set to "auto".

        The 'nticks' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["nticks"]

    @nticks.setter
    def nticks(self, val):
        self["nticks"] = val

    # overlaying
    # ----------
    @property
    def overlaying(self):
        """
        If set a same-letter axis id, this axis is overlaid on top of
        the corresponding same-letter axis, with traces and axes
        visible for both axes. If False, this axis does not overlay any
        same-letter axes. In this case, for axes with overlapping
        domains only the highest-numbered axis will be visible.

        The 'overlaying' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['free']
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$',
                '^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["overlaying"]

    @overlaying.setter
    def overlaying(self, val):
        self["overlaying"] = val

    # position
    # --------
    @property
    def position(self):
        """
        Sets the position of this axis in the plotting space (in
        normalized coordinates). Only has an effect if `anchor` is set
        to "free".

        The 'position' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["position"]

    @position.setter
    def position(self, val):
        self["position"] = val

    # range
    # -----
    @property
    def range(self):
        """
            Sets the range of this axis. If the axis `type` is "log", then
            you must take the log of your desired range (e.g. to set the
            range from 1 to 100, set the range from 0 to 2). If the axis
            `type` is "date", it should be date strings, like date data,
            though Date objects and unix milliseconds will be accepted and
            converted to strings. If the axis `type` is "category", it
            should be numbers, using the scale where each category is
            assigned a serial number from zero in the order it appears.
            Leaving either or both elements `null` impacts the default
            `autorange`.

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

    # rangebreaks
    # -----------
    @property
    def rangebreaks(self):
        """
        The 'rangebreaks' property is a tuple of instances of
        Rangebreak that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.xaxis.Rangebreak
          - A list or tuple of dicts of string/value properties that
            will be passed to the Rangebreak constructor

            Supported dict properties:

                bounds
                    Sets the lower and upper bounds of this axis
                    rangebreak. Can be used with `pattern`.
                dvalue
                    Sets the size of each `values` item. The
                    default is one day in milliseconds.
                enabled
                    Determines whether this axis rangebreak is
                    enabled or disabled. Please note that
                    `rangebreaks` only work for "date" axis type.
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
                pattern
                    Determines a pattern on the time line that
                    generates breaks. If *day of week* - days of
                    the week in English e.g. 'Sunday' or `sun`
                    (matching is case-insensitive and considers
                    only the first three characters), as well as
                    Sunday-based integers between 0 and 6. If
                    "hour" - hour (24-hour clock) as decimal
                    numbers between 0 and 24. for more info.
                    Examples: - { pattern: 'day of week', bounds:
                    [6, 1] }  or simply { bounds: ['sat', 'mon'] }
                    breaks from Saturday to Monday (i.e. skips the
                    weekends). - { pattern: 'hour', bounds: [17, 8]
                    }   breaks from 5pm to 8am (i.e. skips non-work
                    hours).
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
                values
                    Sets the coordinate values corresponding to the
                    rangebreaks. An alternative to `bounds`. Use
                    `dvalue` to set the size of the values along
                    the axis.

        Returns
        -------
        tuple[plotly.graph_objs.layout.xaxis.Rangebreak]
        """
        return self["rangebreaks"]

    @rangebreaks.setter
    def rangebreaks(self, val):
        self["rangebreaks"] = val

    # rangebreakdefaults
    # ------------------
    @property
    def rangebreakdefaults(self):
        """
        When used in a template (as
        layout.template.layout.xaxis.rangebreakdefaults), sets the
        default property values to use for elements of
        layout.xaxis.rangebreaks

        The 'rangebreakdefaults' property is an instance of Rangebreak
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Rangebreak`
          - A dict of string/value properties that will be passed
            to the Rangebreak constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Rangebreak
        """
        return self["rangebreakdefaults"]

    @rangebreakdefaults.setter
    def rangebreakdefaults(self, val):
        self["rangebreakdefaults"] = val

    # rangemode
    # ---------
    @property
    def rangemode(self):
        """
        If "normal", the range is computed in relation to the extrema
        of the input data. If *tozero*`, the range extends to 0,
        regardless of the input data If "nonnegative", the range is
        non-negative, regardless of the input data. Applies only to
        linear axes.

        The 'rangemode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['normal', 'tozero', 'nonnegative']

        Returns
        -------
        Any
        """
        return self["rangemode"]

    @rangemode.setter
    def rangemode(self, val):
        self["rangemode"] = val

    # rangeselector
    # -------------
    @property
    def rangeselector(self):
        """
        The 'rangeselector' property is an instance of Rangeselector
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Rangeselector`
          - A dict of string/value properties that will be passed
            to the Rangeselector constructor

            Supported dict properties:

                activecolor
                    Sets the background color of the active range
                    selector button.
                bgcolor
                    Sets the background color of the range selector
                    buttons.
                bordercolor
                    Sets the color of the border enclosing the
                    range selector.
                borderwidth
                    Sets the width (in px) of the border enclosing
                    the range selector.
                buttons
                    Sets the specifications for each buttons. By
                    default, a range selector comes with no
                    buttons.
                buttondefaults
                    When used in a template (as layout.template.lay
                    out.xaxis.rangeselector.buttondefaults), sets
                    the default property values to use for elements
                    of layout.xaxis.rangeselector.buttons
                font
                    Sets the font of the range selector button
                    text.
                visible
                    Determines whether or not this range selector
                    is visible. Note that range selectors are only
                    available for x axes of `type` set to or auto-
                    typed to "date".
                x
                    Sets the x position (in normalized coordinates)
                    of the range selector.
                xanchor
                    Sets the range selector's horizontal position
                    anchor. This anchor binds the `x` position to
                    the "left", "center" or "right" of the range
                    selector.
                y
                    Sets the y position (in normalized coordinates)
                    of the range selector.
                yanchor
                    Sets the range selector's vertical position
                    anchor This anchor binds the `y` position to
                    the "top", "middle" or "bottom" of the range
                    selector.

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Rangeselector
        """
        return self["rangeselector"]

    @rangeselector.setter
    def rangeselector(self, val):
        self["rangeselector"] = val

    # rangeslider
    # -----------
    @property
    def rangeslider(self):
        """
        The 'rangeslider' property is an instance of Rangeslider
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Rangeslider`
          - A dict of string/value properties that will be passed
            to the Rangeslider constructor

            Supported dict properties:

                autorange
                    Determines whether or not the range slider
                    range is computed in relation to the input
                    data. If `range` is provided, then `autorange`
                    is set to False.
                bgcolor
                    Sets the background color of the range slider.
                bordercolor
                    Sets the border color of the range slider.
                borderwidth
                    Sets the border width of the range slider.
                range
                    Sets the range of the range slider. If not set,
                    defaults to the full xaxis range. If the axis
                    `type` is "log", then you must take the log of
                    your desired range. If the axis `type` is
                    "date", it should be date strings, like date
                    data, though Date objects and unix milliseconds
                    will be accepted and converted to strings. If
                    the axis `type` is "category", it should be
                    numbers, using the scale where each category is
                    assigned a serial number from zero in the order
                    it appears.
                thickness
                    The height of the range slider as a fraction of
                    the total plot area height.
                visible
                    Determines whether or not the range slider will
                    be visible. If visible, perpendicular axes will
                    be set to `fixedrange`
                yaxis
                    :class:`plotly.graph_objects.layout.xaxis.range
                    slider.YAxis` instance or dict with compatible
                    properties

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Rangeslider
        """
        return self["rangeslider"]

    @rangeslider.setter
    def rangeslider(self, val):
        self["rangeslider"] = val

    # scaleanchor
    # -----------
    @property
    def scaleanchor(self):
        """
        If set to another axis id (e.g. `x2`, `y`), the range of this
        axis changes together with the range of the corresponding axis
        such that the scale of pixels per unit is in a constant ratio.
        Both axes are still zoomable, but when you zoom one, the other
        will zoom the same amount, keeping a fixed midpoint.
        `constrain` and `constraintoward` determine how we enforce the
        constraint. You can chain these, ie `yaxis: {scaleanchor: *x*},
        xaxis2: {scaleanchor: *y*}` but you can only link axes of the
        same `type`. The linked axis can have the opposite letter (to
        constrain the aspect ratio) or the same letter (to match scales
        across subplots). Loops (`yaxis: {scaleanchor: *x*}, xaxis:
        {scaleanchor: *y*}` or longer) are redundant and the last
        constraint encountered will be ignored to avoid possible
        inconsistent constraints via `scaleratio`. Note that setting
        axes simultaneously in both a `scaleanchor` and a `matches`
        constraint is currently forbidden. Setting `false` allows to
        remove a default constraint (occasionally, you may need to
        prevent a default `scaleanchor` constraint from being applied,
        eg. when having an image trace `yaxis: {scaleanchor: "x"}` is
        set automatically in order for pixels to be rendered as
        squares, setting `yaxis: {scaleanchor: false}` allows to remove
        the constraint).

        The 'scaleanchor' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [False]
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$',
                '^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["scaleanchor"]

    @scaleanchor.setter
    def scaleanchor(self, val):
        self["scaleanchor"] = val

    # scaleratio
    # ----------
    @property
    def scaleratio(self):
        """
        If this axis is linked to another by `scaleanchor`, this
        determines the pixel to unit scale ratio. For example, if this
        value is 10, then every unit on this axis spans 10 times the
        number of pixels as a unit on the linked axis. Use this for
        example to create an elevation profile where the vertical scale
        is exaggerated a fixed amount with respect to the horizontal.

        The 'scaleratio' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["scaleratio"]

    @scaleratio.setter
    def scaleratio(self, val):
        self["scaleratio"] = val

    # separatethousands
    # -----------------
    @property
    def separatethousands(self):
        """
        If "true", even 4-digit integers are separated

        The 'separatethousands' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["separatethousands"]

    @separatethousands.setter
    def separatethousands(self, val):
        self["separatethousands"] = val

    # showdividers
    # ------------
    @property
    def showdividers(self):
        """
        Determines whether or not a dividers are drawn between the
        category levels of this axis. Only has an effect on
        "multicategory" axes.

        The 'showdividers' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showdividers"]

    @showdividers.setter
    def showdividers(self, val):
        self["showdividers"] = val

    # showexponent
    # ------------
    @property
    def showexponent(self):
        """
        If "all", all exponents are shown besides their significands.
        If "first", only the exponent of the first tick is shown. If
        "last", only the exponent of the last tick is shown. If "none",
        no exponents appear.

        The 'showexponent' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['all', 'first', 'last', 'none']

        Returns
        -------
        Any
        """
        return self["showexponent"]

    @showexponent.setter
    def showexponent(self, val):
        self["showexponent"] = val

    # showgrid
    # --------
    @property
    def showgrid(self):
        """
        Determines whether or not grid lines are drawn. If True, the
        grid lines are drawn at every tick mark.

        The 'showgrid' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showgrid"]

    @showgrid.setter
    def showgrid(self, val):
        self["showgrid"] = val

    # showline
    # --------
    @property
    def showline(self):
        """
        Determines whether or not a line bounding this axis is drawn.

        The 'showline' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showline"]

    @showline.setter
    def showline(self, val):
        self["showline"] = val

    # showspikes
    # ----------
    @property
    def showspikes(self):
        """
        Determines whether or not spikes (aka droplines) are drawn for
        this axis. Note: This only takes affect when hovermode =
        closest

        The 'showspikes' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showspikes"]

    @showspikes.setter
    def showspikes(self, val):
        self["showspikes"] = val

    # showticklabels
    # --------------
    @property
    def showticklabels(self):
        """
        Determines whether or not the tick labels are drawn.

        The 'showticklabels' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showticklabels"]

    @showticklabels.setter
    def showticklabels(self, val):
        self["showticklabels"] = val

    # showtickprefix
    # --------------
    @property
    def showtickprefix(self):
        """
        If "all", all tick labels are displayed with a prefix. If
        "first", only the first tick is displayed with a prefix. If
        "last", only the last tick is displayed with a suffix. If
        "none", tick prefixes are hidden.

        The 'showtickprefix' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['all', 'first', 'last', 'none']

        Returns
        -------
        Any
        """
        return self["showtickprefix"]

    @showtickprefix.setter
    def showtickprefix(self, val):
        self["showtickprefix"] = val

    # showticksuffix
    # --------------
    @property
    def showticksuffix(self):
        """
        Same as `showtickprefix` but for tick suffixes.

        The 'showticksuffix' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['all', 'first', 'last', 'none']

        Returns
        -------
        Any
        """
        return self["showticksuffix"]

    @showticksuffix.setter
    def showticksuffix(self, val):
        self["showticksuffix"] = val

    # side
    # ----
    @property
    def side(self):
        """
        Determines whether a x (y) axis is positioned at the "bottom"
        ("left") or "top" ("right") of the plotting area.

        The 'side' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'bottom', 'left', 'right']

        Returns
        -------
        Any
        """
        return self["side"]

    @side.setter
    def side(self, val):
        self["side"] = val

    # spikecolor
    # ----------
    @property
    def spikecolor(self):
        """
        Sets the spike color. If undefined, will use the series color

        The 'spikecolor' property is a color and may be specified as:
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
        return self["spikecolor"]

    @spikecolor.setter
    def spikecolor(self, val):
        self["spikecolor"] = val

    # spikedash
    # ---------
    @property
    def spikedash(self):
        """
        Sets the dash style of lines. Set to a dash type string
        ("solid", "dot", "dash", "longdash", "dashdot", or
        "longdashdot") or a dash length list in px (eg
        "5px,10px,2px,2px").

        The 'spikedash' property is an enumeration that may be specified as:
          - One of the following dash styles:
                ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
          - A string containing a dash length list in pixels or percentages
                (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)

        Returns
        -------
        str
        """
        return self["spikedash"]

    @spikedash.setter
    def spikedash(self, val):
        self["spikedash"] = val

    # spikemode
    # ---------
    @property
    def spikemode(self):
        """
        Determines the drawing mode for the spike line If "toaxis", the
        line is drawn from the data point to the axis the  series is
        plotted on. If "across", the line is drawn across the entire
        plot area, and supercedes "toaxis". If "marker", then a marker
        dot is drawn on the axis the series is plotted on

        The 'spikemode' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['toaxis', 'across', 'marker'] joined with '+' characters
            (e.g. 'toaxis+across')

        Returns
        -------
        Any
        """
        return self["spikemode"]

    @spikemode.setter
    def spikemode(self, val):
        self["spikemode"] = val

    # spikesnap
    # ---------
    @property
    def spikesnap(self):
        """
        Determines whether spikelines are stuck to the cursor or to the
        closest datapoints.

        The 'spikesnap' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['data', 'cursor', 'hovered data']

        Returns
        -------
        Any
        """
        return self["spikesnap"]

    @spikesnap.setter
    def spikesnap(self, val):
        self["spikesnap"] = val

    # spikethickness
    # --------------
    @property
    def spikethickness(self):
        """
        Sets the width (in px) of the zero line.

        The 'spikethickness' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["spikethickness"]

    @spikethickness.setter
    def spikethickness(self, val):
        self["spikethickness"] = val

    # tick0
    # -----
    @property
    def tick0(self):
        """
        Sets the placement of the first tick on this axis. Use with
        `dtick`. If the axis `type` is "log", then you must take the
        log of your starting tick (e.g. to set the starting tick to
        100, set the `tick0` to 2) except when `dtick`=*L<f>* (see
        `dtick` for more info). If the axis `type` is "date", it should
        be a date string, like date data. If the axis `type` is
        "category", it should be a number, using the scale where each
        category is assigned a serial number from zero in the order it
        appears.

        The 'tick0' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["tick0"]

    @tick0.setter
    def tick0(self, val):
        self["tick0"] = val

    # tickangle
    # ---------
    @property
    def tickangle(self):
        """
        Sets the angle of the tick labels with respect to the
        horizontal. For example, a `tickangle` of -90 draws the tick
        labels vertically.

        The 'tickangle' property is a angle (in degrees) that may be
        specified as a number between -180 and 180.
        Numeric values outside this range are converted to the equivalent value
        (e.g. 270 is converted to -90).

        Returns
        -------
        int|float
        """
        return self["tickangle"]

    @tickangle.setter
    def tickangle(self, val):
        self["tickangle"] = val

    # tickcolor
    # ---------
    @property
    def tickcolor(self):
        """
        Sets the tick color.

        The 'tickcolor' property is a color and may be specified as:
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
        return self["tickcolor"]

    @tickcolor.setter
    def tickcolor(self, val):
        self["tickcolor"] = val

    # tickfont
    # --------
    @property
    def tickfont(self):
        """
        Sets the tick font.

        The 'tickfont' property is an instance of Tickfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Tickfont`
          - A dict of string/value properties that will be passed
            to the Tickfont constructor

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
        plotly.graph_objs.layout.xaxis.Tickfont
        """
        return self["tickfont"]

    @tickfont.setter
    def tickfont(self, val):
        self["tickfont"] = val

    # tickformat
    # ----------
    @property
    def tickformat(self):
        """
        Sets the tick label formatting rule using d3 formatting mini-
        languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format. And for
        dates see: https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format. We add two items to d3's date
        formatter: "%h" for half of the year as a decimal number as
        well as "%{n}f" for fractional seconds with n digits. For
        example, *2016-10-13 09:15:23.456* with tickformat
        "%H~%M~%S.%2f" would display "09~15~23.46"

        The 'tickformat' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["tickformat"]

    @tickformat.setter
    def tickformat(self, val):
        self["tickformat"] = val

    # tickformatstops
    # ---------------
    @property
    def tickformatstops(self):
        """
        The 'tickformatstops' property is a tuple of instances of
        Tickformatstop that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.xaxis.Tickformatstop
          - A list or tuple of dicts of string/value properties that
            will be passed to the Tickformatstop constructor

            Supported dict properties:

                dtickrange
                    range [*min*, *max*], where "min", "max" -
                    dtick values which describe some zoom level, it
                    is possible to omit "min" or "max" value by
                    passing "null"
                enabled
                    Determines whether or not this stop is used. If
                    `false`, this stop is ignored even within its
                    `dtickrange`.
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
                value
                    string - dtickformat for described zoom level,
                    the same as "tickformat"

        Returns
        -------
        tuple[plotly.graph_objs.layout.xaxis.Tickformatstop]
        """
        return self["tickformatstops"]

    @tickformatstops.setter
    def tickformatstops(self, val):
        self["tickformatstops"] = val

    # tickformatstopdefaults
    # ----------------------
    @property
    def tickformatstopdefaults(self):
        """
        When used in a template (as
        layout.template.layout.xaxis.tickformatstopdefaults), sets the
        default property values to use for elements of
        layout.xaxis.tickformatstops

        The 'tickformatstopdefaults' property is an instance of Tickformatstop
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Tickformatstop`
          - A dict of string/value properties that will be passed
            to the Tickformatstop constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Tickformatstop
        """
        return self["tickformatstopdefaults"]

    @tickformatstopdefaults.setter
    def tickformatstopdefaults(self, val):
        self["tickformatstopdefaults"] = val

    # ticklabelindex
    # --------------
    @property
    def ticklabelindex(self):
        """
        Only for axes with `type` "date" or "linear". Instead of
        drawing the major tick label, draw the label for the minor tick
        that is n positions away from the major tick. E.g. to always
        draw the label for the minor tick before each major tick,
        choose `ticklabelindex` -1. This is useful for date axes with
        `ticklabelmode` "period" if you want to label the period that
        ends with each major tick instead of the period that begins
        there.

        The 'ticklabelindex' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|numpy.ndarray
        """
        return self["ticklabelindex"]

    @ticklabelindex.setter
    def ticklabelindex(self, val):
        self["ticklabelindex"] = val

    # ticklabelindexsrc
    # -----------------
    @property
    def ticklabelindexsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `ticklabelindex`.

        The 'ticklabelindexsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["ticklabelindexsrc"]

    @ticklabelindexsrc.setter
    def ticklabelindexsrc(self, val):
        self["ticklabelindexsrc"] = val

    # ticklabelmode
    # -------------
    @property
    def ticklabelmode(self):
        """
        Determines where tick labels are drawn with respect to their
        corresponding ticks and grid lines. Only has an effect for axes
        of `type` "date" When set to "period", tick labels are drawn in
        the middle of the period between ticks.

        The 'ticklabelmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['instant', 'period']

        Returns
        -------
        Any
        """
        return self["ticklabelmode"]

    @ticklabelmode.setter
    def ticklabelmode(self, val):
        self["ticklabelmode"] = val

    # ticklabeloverflow
    # -----------------
    @property
    def ticklabeloverflow(self):
        """
        Determines how we handle tick labels that would overflow either
        the graph div or the domain of the axis. The default value for
        inside tick labels is *hide past domain*. Otherwise on
        "category" and "multicategory" axes the default is "allow". In
        other cases the default is *hide past div*.

        The 'ticklabeloverflow' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['allow', 'hide past div', 'hide past domain']

        Returns
        -------
        Any
        """
        return self["ticklabeloverflow"]

    @ticklabeloverflow.setter
    def ticklabeloverflow(self, val):
        self["ticklabeloverflow"] = val

    # ticklabelposition
    # -----------------
    @property
    def ticklabelposition(self):
        """
        Determines where tick labels are drawn with respect to the axis
        Please note that top or bottom has no effect on x axes or when
        `ticklabelmode` is set to "period". Similarly left or right has
        no effect on y axes or when `ticklabelmode` is set to "period".
        Has no effect on "multicategory" axes or when `tickson` is set
        to "boundaries". When used on axes linked by `matches` or
        `scaleanchor`, no extra padding for inside labels would be
        added by autorange, so that the scales could match.

        The 'ticklabelposition' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['outside', 'inside', 'outside top', 'inside top',
                'outside left', 'inside left', 'outside right', 'inside
                right', 'outside bottom', 'inside bottom']

        Returns
        -------
        Any
        """
        return self["ticklabelposition"]

    @ticklabelposition.setter
    def ticklabelposition(self, val):
        self["ticklabelposition"] = val

    # ticklabelshift
    # --------------
    @property
    def ticklabelshift(self):
        """
        Shifts the tick labels by the specified number of pixels in
        parallel to the axis. Positive values move the labels in the
        positive direction of the axis.

        The 'ticklabelshift' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)

        Returns
        -------
        int
        """
        return self["ticklabelshift"]

    @ticklabelshift.setter
    def ticklabelshift(self, val):
        self["ticklabelshift"] = val

    # ticklabelstandoff
    # -----------------
    @property
    def ticklabelstandoff(self):
        """
        Sets the standoff distance (in px) between the axis tick labels
        and their default position. A positive `ticklabelstandoff`
        moves the labels farther away from the plot area if
        `ticklabelposition` is "outside", and deeper into the plot area
        if `ticklabelposition` is "inside". A negative
        `ticklabelstandoff` works in the opposite direction, moving
        outside ticks towards the plot area and inside ticks towards
        the outside. If the negative value is large enough, inside
        ticks can even end up outside and vice versa.

        The 'ticklabelstandoff' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)

        Returns
        -------
        int
        """
        return self["ticklabelstandoff"]

    @ticklabelstandoff.setter
    def ticklabelstandoff(self, val):
        self["ticklabelstandoff"] = val

    # ticklabelstep
    # -------------
    @property
    def ticklabelstep(self):
        """
        Sets the spacing between tick labels as compared to the spacing
        between ticks. A value of 1 (default) means each tick gets a
        label. A value of 2 means shows every 2nd label. A larger value
        n means only every nth tick is labeled. `tick0` determines
        which labels are shown. Not implemented for axes with `type`
        "log" or "multicategory", or when `tickmode` is "array".

        The 'ticklabelstep' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["ticklabelstep"]

    @ticklabelstep.setter
    def ticklabelstep(self, val):
        self["ticklabelstep"] = val

    # ticklen
    # -------
    @property
    def ticklen(self):
        """
        Sets the tick length (in px).

        The 'ticklen' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["ticklen"]

    @ticklen.setter
    def ticklen(self, val):
        self["ticklen"] = val

    # tickmode
    # --------
    @property
    def tickmode(self):
        """
        Sets the tick mode for this axis. If "auto", the number of
        ticks is set via `nticks`. If "linear", the placement of the
        ticks is determined by a starting position `tick0` and a tick
        step `dtick` ("linear" is the default value if `tick0` and
        `dtick` are provided). If "array", the placement of the ticks
        is set via `tickvals` and the tick text is `ticktext`. ("array"
        is the default value if `tickvals` is provided). If "sync", the
        number of ticks will sync with the overlayed axis set by
        `overlaying` property.

        The 'tickmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['auto', 'linear', 'array', 'sync']

        Returns
        -------
        Any
        """
        return self["tickmode"]

    @tickmode.setter
    def tickmode(self, val):
        self["tickmode"] = val

    # tickprefix
    # ----------
    @property
    def tickprefix(self):
        """
        Sets a tick label prefix.

        The 'tickprefix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["tickprefix"]

    @tickprefix.setter
    def tickprefix(self, val):
        self["tickprefix"] = val

    # ticks
    # -----
    @property
    def ticks(self):
        """
        Determines whether ticks are drawn or not. If "", this axis'
        ticks are not drawn. If "outside" ("inside"), this axis' are
        drawn outside (inside) the axis lines.

        The 'ticks' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['outside', 'inside', '']

        Returns
        -------
        Any
        """
        return self["ticks"]

    @ticks.setter
    def ticks(self, val):
        self["ticks"] = val

    # tickson
    # -------
    @property
    def tickson(self):
        """
        Determines where ticks and grid lines are drawn with respect to
        their corresponding tick labels. Only has an effect for axes of
        `type` "category" or "multicategory". When set to "boundaries",
        ticks and grid lines are drawn half a category to the
        left/bottom of labels.

        The 'tickson' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['labels', 'boundaries']

        Returns
        -------
        Any
        """
        return self["tickson"]

    @tickson.setter
    def tickson(self, val):
        self["tickson"] = val

    # ticksuffix
    # ----------
    @property
    def ticksuffix(self):
        """
        Sets a tick label suffix.

        The 'ticksuffix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["ticksuffix"]

    @ticksuffix.setter
    def ticksuffix(self, val):
        self["ticksuffix"] = val

    # ticktext
    # --------
    @property
    def ticktext(self):
        """
        Sets the text displayed at the ticks position via `tickvals`.
        Only has an effect if `tickmode` is set to "array". Used with
        `tickvals`.

        The 'ticktext' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["ticktext"]

    @ticktext.setter
    def ticktext(self, val):
        self["ticktext"] = val

    # ticktextsrc
    # -----------
    @property
    def ticktextsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `ticktext`.

        The 'ticktextsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["ticktextsrc"]

    @ticktextsrc.setter
    def ticktextsrc(self, val):
        self["ticktextsrc"] = val

    # tickvals
    # --------
    @property
    def tickvals(self):
        """
        Sets the values at which ticks on this axis appear. Only has an
        effect if `tickmode` is set to "array". Used with `ticktext`.

        The 'tickvals' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["tickvals"]

    @tickvals.setter
    def tickvals(self, val):
        self["tickvals"] = val

    # tickvalssrc
    # -----------
    @property
    def tickvalssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `tickvals`.

        The 'tickvalssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["tickvalssrc"]

    @tickvalssrc.setter
    def tickvalssrc(self, val):
        self["tickvalssrc"] = val

    # tickwidth
    # ---------
    @property
    def tickwidth(self):
        """
        Sets the tick width (in px).

        The 'tickwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["tickwidth"]

    @tickwidth.setter
    def tickwidth(self, val):
        self["tickwidth"] = val

    # title
    # -----
    @property
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

            Supported dict properties:

                font
                    Sets this axis' title font. Note that the
                    title's font used to be customized by the now
                    deprecated `titlefont` attribute.
                standoff
                    Sets the standoff distance (in px) between the
                    axis labels and the title text The default
                    value is a function of the axis tick labels,
                    the title `font.size` and the axis `linewidth`.
                    Note that the axis title position is always
                    constrained within the margins, so the actual
                    standoff distance is always less than the set
                    or default value. By setting `standoff` and
                    turning on `automargin`, plotly.js will push
                    the margins to fit the axis title at given
                    standoff distance.
                text
                    Sets the title of this axis. Note that before
                    the existence of `title.text`, the title's
                    contents used to be defined as the `title`
                    attribute itself. This behavior has been
                    deprecated.

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Title
        """
        return self["title"]

    @title.setter
    def title(self, val):
        self["title"] = val

    # titlefont
    # ---------
    @property
    def titlefont(self):
        """
        Deprecated: Please use layout.xaxis.title.font instead. Sets
        this axis' title font. Note that the title's font used to be
        customized by the now deprecated `titlefont` attribute.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.title.Font`
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

        """
        return self["titlefont"]

    @titlefont.setter
    def titlefont(self, val):
        self["titlefont"] = val

    # type
    # ----
    @property
    def type(self):
        """
        Sets the axis type. By default, plotly attempts to determined
        the axis type by looking into the data of the traces that
        referenced the axis in question.

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['-', 'linear', 'log', 'date', 'category',
                'multicategory']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    # uirevision
    # ----------
    @property
    def uirevision(self):
        """
        Controls persistence of user-driven changes in axis `range`,
        `autorange`, and `title` if in `editable: true` configuration.
        Defaults to `layout.uirevision`.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        A single toggle to hide the axis while preserving interaction
        like dragging. Default is true when a cheater plot is present
        on the axis, otherwise false

        The 'visible' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    # zeroline
    # --------
    @property
    def zeroline(self):
        """
        Determines whether or not a line is drawn at along the 0 value
        of this axis. If True, the zero line is drawn on top of the
        grid lines.

        The 'zeroline' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["zeroline"]

    @zeroline.setter
    def zeroline(self, val):
        self["zeroline"] = val

    # zerolinecolor
    # -------------
    @property
    def zerolinecolor(self):
        """
        Sets the line color of the zero line.

        The 'zerolinecolor' property is a color and may be specified as:
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
        return self["zerolinecolor"]

    @zerolinecolor.setter
    def zerolinecolor(self, val):
        self["zerolinecolor"] = val

    # zerolinewidth
    # -------------
    @property
    def zerolinewidth(self):
        """
        Sets the width (in px) of the zero line.

        The 'zerolinewidth' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["zerolinewidth"]

    @zerolinewidth.setter
    def zerolinewidth(self, val):
        self["zerolinewidth"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        anchor
            If set to an opposite-letter axis id (e.g. `x2`, `y`),
            this axis is bound to the corresponding opposite-letter
            axis. If set to "free", this axis' position is
            determined by `position`.
        automargin
            Determines whether long tick labels automatically grow
            the figure margins.
        autorange
            Determines whether or not the range of this axis is
            computed in relation to the input data. See `rangemode`
            for more info. If `range` is provided and it has a
            value for both the lower and upper bound, `autorange`
            is set to False. Using "min" applies autorange only to
            set the minimum. Using "max" applies autorange only to
            set the maximum. Using *min reversed* applies autorange
            only to set the minimum on a reversed axis. Using *max
            reversed* applies autorange only to set the maximum on
            a reversed axis. Using "reversed" applies autorange on
            both ends and reverses the axis direction.
        autorangeoptions
            :class:`plotly.graph_objects.layout.xaxis.Autorangeopti
            ons` instance or dict with compatible properties
        autotickangles
            When `tickangle` is set to "auto", it will be set to
            the first angle in this array that is large enough to
            prevent label overlap.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. Defaults to
            layout.autotypenumbers.
        calendar
            Sets the calendar system to use for `range` and `tick0`
            if this is a date axis. This does not set the calendar
            for interpreting data on this axis, that's specified in
            the trace or via the global `layout.calendar`
        categoryarray
            Sets the order in which categories on this axis appear.
            Only has an effect if `categoryorder` is set to
            "array". Used with `categoryorder`.
        categoryarraysrc
            Sets the source reference on Chart Studio Cloud for
            `categoryarray`.
        categoryorder
            Specifies the ordering logic for the case of
            categorical variables. By default, plotly uses "trace",
            which specifies the order that is present in the data
            supplied. Set `categoryorder` to *category ascending*
            or *category descending* if order should be determined
            by the alphanumerical order of the category names. Set
            `categoryorder` to "array" to derive the ordering from
            the attribute `categoryarray`. If a category is not
            found in the `categoryarray` array, the sorting
            behavior for that attribute will be identical to the
            "trace" mode. The unspecified categories will follow
            the categories in `categoryarray`. Set `categoryorder`
            to *total ascending* or *total descending* if order
            should be determined by the numerical order of the
            values. Similarly, the order can be determined by the
            min, max, sum, mean, geometric mean or median of all
            the values.
        color
            Sets default for all colors associated with this axis
            all at once: line, font, tick, and grid colors. Grid
            color is lightened by blending this with the plot
            background Individual pieces can override this.
        constrain
            If this axis needs to be compressed (either due to its
            own `scaleanchor` and `scaleratio` or those of the
            other axis), determines how that happens: by increasing
            the "range", or by decreasing the "domain". Default is
            "domain" for axes containing image traces, "range"
            otherwise.
        constraintoward
            If this axis needs to be compressed (either due to its
            own `scaleanchor` and `scaleratio` or those of the
            other axis), determines which direction we push the
            originally specified plot area. Options are "left",
            "center" (default), and "right" for x axes, and "top",
            "middle" (default), and "bottom" for y axes.
        dividercolor
            Sets the color of the dividers Only has an effect on
            "multicategory" axes.
        dividerwidth
            Sets the width (in px) of the dividers Only has an
            effect on "multicategory" axes.
        domain
            Sets the domain of this axis (in plot fraction).
        dtick
            Sets the step in-between ticks on this axis. Use with
            `tick0`. Must be a positive number, or special strings
            available to "log" and "date" axes. If the axis `type`
            is "log", then ticks are set every 10^(n*dtick) where n
            is the tick number. For example, to set a tick mark at
            1, 10, 100, 1000, ... set dtick to 1. To set tick marks
            at 1, 100, 10000, ... set dtick to 2. To set tick marks
            at 1, 5, 25, 125, 625, 3125, ... set dtick to
            log_10(5), or 0.69897000433. "log" has several special
            values; "L<f>", where `f` is a positive number, gives
            ticks linearly spaced in value (but not position). For
            example `tick0` = 0.1, `dtick` = "L0.5" will put ticks
            at 0.1, 0.6, 1.1, 1.6 etc. To show powers of 10 plus
            small digits between, use "D1" (all digits) or "D2"
            (only 2 and 5). `tick0` is ignored for "D1" and "D2".
            If the axis `type` is "date", then you must convert the
            time to milliseconds. For example, to set the interval
            between ticks to one day, set `dtick` to 86400000.0.
            "date" also has special values "M<n>" gives ticks
            spaced by a number of months. `n` must be a positive
            integer. To set ticks on the 15th of every third month,
            set `tick0` to "2000-01-15" and `dtick` to "M3". To set
            ticks every 4 years, set `dtick` to "M48"
        exponentformat
            Determines a formatting rule for the tick exponents.
            For example, consider the number 1,000,000,000. If
            "none", it appears as 1,000,000,000. If "e", 1e+9. If
            "E", 1E+9. If "power", 1x10^9 (with 9 in a super
            script). If "SI", 1G. If "B", 1B.
        fixedrange
            Determines whether or not this axis is zoom-able. If
            true, then zoom is disabled.
        gridcolor
            Sets the color of the grid lines.
        griddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        gridwidth
            Sets the width (in px) of the grid lines.
        hoverformat
            Sets the hover text formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display "09~15~23.46"
        insiderange
            Could be used to set the desired inside range of this
            axis (excluding the labels) when `ticklabelposition` of
            the anchored axis has "inside". Not implemented for
            axes with `type` "log". This would be ignored when
            `range` is provided.
        labelalias
            Replacement text for specific tick or hover labels. For
            example using {US: 'USA', CA: 'Canada'} changes US to
            USA and CA to Canada. The labels we would have shown
            must match the keys exactly, after adding any
            tickprefix or ticksuffix. For negative numbers the
            minus sign symbol used (U+2212) is wider than the
            regular ascii dash. That means you need to use −1
            instead of -1. labelalias can be used with any axis
            type, and both keys (if needed) and values (if desired)
            can include html-like tags or MathJax.
        layer
            Sets the layer on which this axis is displayed. If
            *above traces*, this axis is displayed above all the
            subplot's traces If *below traces*, this axis is
            displayed below all the subplot's traces, but above the
            grid lines. Useful when used together with scatter-like
            traces with `cliponaxis` set to False to show markers
            and/or text nodes above this axis.
        linecolor
            Sets the axis line color.
        linewidth
            Sets the width (in px) of the axis line.
        matches
            If set to another axis id (e.g. `x2`, `y`), the range
            of this axis will match the range of the corresponding
            axis in data-coordinates space. Moreover, matching axes
            share auto-range values, category lists and histogram
            auto-bins. Note that setting axes simultaneously in
            both a `scaleanchor` and a `matches` constraint is
            currently forbidden. Moreover, note that matching axes
            must have the same `type`.
        maxallowed
            Determines the maximum range of this axis.
        minallowed
            Determines the minimum range of this axis.
        minexponent
            Hide SI prefix for 10^n if |n| is below this number.
            This only has an effect when `tickformat` is "SI" or
            "B".
        minor
            :class:`plotly.graph_objects.layout.xaxis.Minor`
            instance or dict with compatible properties
        mirror
            Determines if the axis lines or/and ticks are mirrored
            to the opposite side of the plotting area. If True, the
            axis lines are mirrored. If "ticks", the axis lines and
            ticks are mirrored. If False, mirroring is disable. If
            "all", axis lines are mirrored on all shared-axes
            subplots. If "allticks", axis lines and ticks are
            mirrored on all shared-axes subplots.
        nticks
            Specifies the maximum number of ticks for the
            particular axis. The actual number of ticks will be
            chosen automatically to be less than or equal to
            `nticks`. Has an effect only if `tickmode` is set to
            "auto".
        overlaying
            If set a same-letter axis id, this axis is overlaid on
            top of the corresponding same-letter axis, with traces
            and axes visible for both axes. If False, this axis
            does not overlay any same-letter axes. In this case,
            for axes with overlapping domains only the highest-
            numbered axis will be visible.
        position
            Sets the position of this axis in the plotting space
            (in normalized coordinates). Only has an effect if
            `anchor` is set to "free".
        range
            Sets the range of this axis. If the axis `type` is
            "log", then you must take the log of your desired range
            (e.g. to set the range from 1 to 100, set the range
            from 0 to 2). If the axis `type` is "date", it should
            be date strings, like date data, though Date objects
            and unix milliseconds will be accepted and converted to
            strings. If the axis `type` is "category", it should be
            numbers, using the scale where each category is
            assigned a serial number from zero in the order it
            appears. Leaving either or both elements `null` impacts
            the default `autorange`.
        rangebreaks
            A tuple of
            :class:`plotly.graph_objects.layout.xaxis.Rangebreak`
            instances or dicts with compatible properties
        rangebreakdefaults
            When used in a template (as
            layout.template.layout.xaxis.rangebreakdefaults), sets
            the default property values to use for elements of
            layout.xaxis.rangebreaks
        rangemode
            If "normal", the range is computed in relation to the
            extrema of the input data. If *tozero*`, the range
            extends to 0, regardless of the input data If
            "nonnegative", the range is non-negative, regardless of
            the input data. Applies only to linear axes.
        rangeselector
            :class:`plotly.graph_objects.layout.xaxis.Rangeselector
            ` instance or dict with compatible properties
        rangeslider
            :class:`plotly.graph_objects.layout.xaxis.Rangeslider`
            instance or dict with compatible properties
        scaleanchor
            If set to another axis id (e.g. `x2`, `y`), the range
            of this axis changes together with the range of the
            corresponding axis such that the scale of pixels per
            unit is in a constant ratio. Both axes are still
            zoomable, but when you zoom one, the other will zoom
            the same amount, keeping a fixed midpoint. `constrain`
            and `constraintoward` determine how we enforce the
            constraint. You can chain these, ie `yaxis:
            {scaleanchor: *x*}, xaxis2: {scaleanchor: *y*}` but you
            can only link axes of the same `type`. The linked axis
            can have the opposite letter (to constrain the aspect
            ratio) or the same letter (to match scales across
            subplots). Loops (`yaxis: {scaleanchor: *x*}, xaxis:
            {scaleanchor: *y*}` or longer) are redundant and the
            last constraint encountered will be ignored to avoid
            possible inconsistent constraints via `scaleratio`.
            Note that setting axes simultaneously in both a
            `scaleanchor` and a `matches` constraint is currently
            forbidden. Setting `false` allows to remove a default
            constraint (occasionally, you may need to prevent a
            default `scaleanchor` constraint from being applied,
            eg. when having an image trace `yaxis: {scaleanchor:
            "x"}` is set automatically in order for pixels to be
            rendered as squares, setting `yaxis: {scaleanchor:
            false}` allows to remove the constraint).
        scaleratio
            If this axis is linked to another by `scaleanchor`,
            this determines the pixel to unit scale ratio. For
            example, if this value is 10, then every unit on this
            axis spans 10 times the number of pixels as a unit on
            the linked axis. Use this for example to create an
            elevation profile where the vertical scale is
            exaggerated a fixed amount with respect to the
            horizontal.
        separatethousands
            If "true", even 4-digit integers are separated
        showdividers
            Determines whether or not a dividers are drawn between
            the category levels of this axis. Only has an effect on
            "multicategory" axes.
        showexponent
            If "all", all exponents are shown besides their
            significands. If "first", only the exponent of the
            first tick is shown. If "last", only the exponent of
            the last tick is shown. If "none", no exponents appear.
        showgrid
            Determines whether or not grid lines are drawn. If
            True, the grid lines are drawn at every tick mark.
        showline
            Determines whether or not a line bounding this axis is
            drawn.
        showspikes
            Determines whether or not spikes (aka droplines) are
            drawn for this axis. Note: This only takes affect when
            hovermode = closest
        showticklabels
            Determines whether or not the tick labels are drawn.
        showtickprefix
            If "all", all tick labels are displayed with a prefix.
            If "first", only the first tick is displayed with a
            prefix. If "last", only the last tick is displayed with
            a suffix. If "none", tick prefixes are hidden.
        showticksuffix
            Same as `showtickprefix` but for tick suffixes.
        side
            Determines whether a x (y) axis is positioned at the
            "bottom" ("left") or "top" ("right") of the plotting
            area.
        spikecolor
            Sets the spike color. If undefined, will use the series
            color
        spikedash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        spikemode
            Determines the drawing mode for the spike line If
            "toaxis", the line is drawn from the data point to the
            axis the  series is plotted on. If "across", the line
            is drawn across the entire plot area, and supercedes
            "toaxis". If "marker", then a marker dot is drawn on
            the axis the series is plotted on
        spikesnap
            Determines whether spikelines are stuck to the cursor
            or to the closest datapoints.
        spikethickness
            Sets the width (in px) of the zero line.
        tick0
            Sets the placement of the first tick on this axis. Use
            with `dtick`. If the axis `type` is "log", then you
            must take the log of your starting tick (e.g. to set
            the starting tick to 100, set the `tick0` to 2) except
            when `dtick`=*L<f>* (see `dtick` for more info). If the
            axis `type` is "date", it should be a date string, like
            date data. If the axis `type` is "category", it should
            be a number, using the scale where each category is
            assigned a serial number from zero in the order it
            appears.
        tickangle
            Sets the angle of the tick labels with respect to the
            horizontal. For example, a `tickangle` of -90 draws the
            tick labels vertically.
        tickcolor
            Sets the tick color.
        tickfont
            Sets the tick font.
        tickformat
            Sets the tick label formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display "09~15~23.46"
        tickformatstops
            A tuple of :class:`plotly.graph_objects.layout.xaxis.Ti
            ckformatstop` instances or dicts with compatible
            properties
        tickformatstopdefaults
            When used in a template (as
            layout.template.layout.xaxis.tickformatstopdefaults),
            sets the default property values to use for elements of
            layout.xaxis.tickformatstops
        ticklabelindex
            Only for axes with `type` "date" or "linear". Instead
            of drawing the major tick label, draw the label for the
            minor tick that is n positions away from the major
            tick. E.g. to always draw the label for the minor tick
            before each major tick, choose `ticklabelindex` -1.
            This is useful for date axes with `ticklabelmode`
            "period" if you want to label the period that ends with
            each major tick instead of the period that begins
            there.
        ticklabelindexsrc
            Sets the source reference on Chart Studio Cloud for
            `ticklabelindex`.
        ticklabelmode
            Determines where tick labels are drawn with respect to
            their corresponding ticks and grid lines. Only has an
            effect for axes of `type` "date" When set to "period",
            tick labels are drawn in the middle of the period
            between ticks.
        ticklabeloverflow
            Determines how we handle tick labels that would
            overflow either the graph div or the domain of the
            axis. The default value for inside tick labels is *hide
            past domain*. Otherwise on "category" and
            "multicategory" axes the default is "allow". In other
            cases the default is *hide past div*.
        ticklabelposition
            Determines where tick labels are drawn with respect to
            the axis Please note that top or bottom has no effect
            on x axes or when `ticklabelmode` is set to "period".
            Similarly left or right has no effect on y axes or when
            `ticklabelmode` is set to "period". Has no effect on
            "multicategory" axes or when `tickson` is set to
            "boundaries". When used on axes linked by `matches` or
            `scaleanchor`, no extra padding for inside labels would
            be added by autorange, so that the scales could match.
        ticklabelshift
            Shifts the tick labels by the specified number of
            pixels in parallel to the axis. Positive values move
            the labels in the positive direction of the axis.
        ticklabelstandoff
            Sets the standoff distance (in px) between the axis
            tick labels and their default position. A positive
            `ticklabelstandoff` moves the labels farther away from
            the plot area if `ticklabelposition` is "outside", and
            deeper into the plot area if `ticklabelposition` is
            "inside". A negative `ticklabelstandoff` works in the
            opposite direction, moving outside ticks towards the
            plot area and inside ticks towards the outside. If the
            negative value is large enough, inside ticks can even
            end up outside and vice versa.
        ticklabelstep
            Sets the spacing between tick labels as compared to the
            spacing between ticks. A value of 1 (default) means
            each tick gets a label. A value of 2 means shows every
            2nd label. A larger value n means only every nth tick
            is labeled. `tick0` determines which labels are shown.
            Not implemented for axes with `type` "log" or
            "multicategory", or when `tickmode` is "array".
        ticklen
            Sets the tick length (in px).
        tickmode
            Sets the tick mode for this axis. If "auto", the number
            of ticks is set via `nticks`. If "linear", the
            placement of the ticks is determined by a starting
            position `tick0` and a tick step `dtick` ("linear" is
            the default value if `tick0` and `dtick` are provided).
            If "array", the placement of the ticks is set via
            `tickvals` and the tick text is `ticktext`. ("array" is
            the default value if `tickvals` is provided). If
            "sync", the number of ticks will sync with the
            overlayed axis set by `overlaying` property.
        tickprefix
            Sets a tick label prefix.
        ticks
            Determines whether ticks are drawn or not. If "", this
            axis' ticks are not drawn. If "outside" ("inside"),
            this axis' are drawn outside (inside) the axis lines.
        tickson
            Determines where ticks and grid lines are drawn with
            respect to their corresponding tick labels. Only has an
            effect for axes of `type` "category" or
            "multicategory". When set to "boundaries", ticks and
            grid lines are drawn half a category to the left/bottom
            of labels.
        ticksuffix
            Sets a tick label suffix.
        ticktext
            Sets the text displayed at the ticks position via
            `tickvals`. Only has an effect if `tickmode` is set to
            "array". Used with `tickvals`.
        ticktextsrc
            Sets the source reference on Chart Studio Cloud for
            `ticktext`.
        tickvals
            Sets the values at which ticks on this axis appear.
            Only has an effect if `tickmode` is set to "array".
            Used with `ticktext`.
        tickvalssrc
            Sets the source reference on Chart Studio Cloud for
            `tickvals`.
        tickwidth
            Sets the tick width (in px).
        title
            :class:`plotly.graph_objects.layout.xaxis.Title`
            instance or dict with compatible properties
        titlefont
            Deprecated: Please use layout.xaxis.title.font instead.
            Sets this axis' title font. Note that the title's font
            used to be customized by the now deprecated `titlefont`
            attribute.
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.
        uirevision
            Controls persistence of user-driven changes in axis
            `range`, `autorange`, and `title` if in `editable:
            true` configuration. Defaults to `layout.uirevision`.
        visible
            A single toggle to hide the axis while preserving
            interaction like dragging. Default is true when a
            cheater plot is present on the axis, otherwise false
        zeroline
            Determines whether or not a line is drawn at along the
            0 value of this axis. If True, the zero line is drawn
            on top of the grid lines.
        zerolinecolor
            Sets the line color of the zero line.
        zerolinewidth
            Sets the width (in px) of the zero line.
        """

    _mapped_properties = {"titlefont": ("title", "font")}

    def __init__(
        self,
        arg=None,
        anchor=None,
        automargin=None,
        autorange=None,
        autorangeoptions=None,
        autotickangles=None,
        autotypenumbers=None,
        calendar=None,
        categoryarray=None,
        categoryarraysrc=None,
        categoryorder=None,
        color=None,
        constrain=None,
        constraintoward=None,
        dividercolor=None,
        dividerwidth=None,
        domain=None,
        dtick=None,
        exponentformat=None,
        fixedrange=None,
        gridcolor=None,
        griddash=None,
        gridwidth=None,
        hoverformat=None,
        insiderange=None,
        labelalias=None,
        layer=None,
        linecolor=None,
        linewidth=None,
        matches=None,
        maxallowed=None,
        minallowed=None,
        minexponent=None,
        minor=None,
        mirror=None,
        nticks=None,
        overlaying=None,
        position=None,
        range=None,
        rangebreaks=None,
        rangebreakdefaults=None,
        rangemode=None,
        rangeselector=None,
        rangeslider=None,
        scaleanchor=None,
        scaleratio=None,
        separatethousands=None,
        showdividers=None,
        showexponent=None,
        showgrid=None,
        showline=None,
        showspikes=None,
        showticklabels=None,
        showtickprefix=None,
        showticksuffix=None,
        side=None,
        spikecolor=None,
        spikedash=None,
        spikemode=None,
        spikesnap=None,
        spikethickness=None,
        tick0=None,
        tickangle=None,
        tickcolor=None,
        tickfont=None,
        tickformat=None,
        tickformatstops=None,
        tickformatstopdefaults=None,
        ticklabelindex=None,
        ticklabelindexsrc=None,
        ticklabelmode=None,
        ticklabeloverflow=None,
        ticklabelposition=None,
        ticklabelshift=None,
        ticklabelstandoff=None,
        ticklabelstep=None,
        ticklen=None,
        tickmode=None,
        tickprefix=None,
        ticks=None,
        tickson=None,
        ticksuffix=None,
        ticktext=None,
        ticktextsrc=None,
        tickvals=None,
        tickvalssrc=None,
        tickwidth=None,
        title=None,
        titlefont=None,
        type=None,
        uirevision=None,
        visible=None,
        zeroline=None,
        zerolinecolor=None,
        zerolinewidth=None,
        **kwargs,
    ):
        """
        Construct a new XAxis object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.XAxis`
        anchor
            If set to an opposite-letter axis id (e.g. `x2`, `y`),
            this axis is bound to the corresponding opposite-letter
            axis. If set to "free", this axis' position is
            determined by `position`.
        automargin
            Determines whether long tick labels automatically grow
            the figure margins.
        autorange
            Determines whether or not the range of this axis is
            computed in relation to the input data. See `rangemode`
            for more info. If `range` is provided and it has a
            value for both the lower and upper bound, `autorange`
            is set to False. Using "min" applies autorange only to
            set the minimum. Using "max" applies autorange only to
            set the maximum. Using *min reversed* applies autorange
            only to set the minimum on a reversed axis. Using *max
            reversed* applies autorange only to set the maximum on
            a reversed axis. Using "reversed" applies autorange on
            both ends and reverses the axis direction.
        autorangeoptions
            :class:`plotly.graph_objects.layout.xaxis.Autorangeopti
            ons` instance or dict with compatible properties
        autotickangles
            When `tickangle` is set to "auto", it will be set to
            the first angle in this array that is large enough to
            prevent label overlap.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. Defaults to
            layout.autotypenumbers.
        calendar
            Sets the calendar system to use for `range` and `tick0`
            if this is a date axis. This does not set the calendar
            for interpreting data on this axis, that's specified in
            the trace or via the global `layout.calendar`
        categoryarray
            Sets the order in which categories on this axis appear.
            Only has an effect if `categoryorder` is set to
            "array". Used with `categoryorder`.
        categoryarraysrc
            Sets the source reference on Chart Studio Cloud for
            `categoryarray`.
        categoryorder
            Specifies the ordering logic for the case of
            categorical variables. By default, plotly uses "trace",
            which specifies the order that is present in the data
            supplied. Set `categoryorder` to *category ascending*
            or *category descending* if order should be determined
            by the alphanumerical order of the category names. Set
            `categoryorder` to "array" to derive the ordering from
            the attribute `categoryarray`. If a category is not
            found in the `categoryarray` array, the sorting
            behavior for that attribute will be identical to the
            "trace" mode. The unspecified categories will follow
            the categories in `categoryarray`. Set `categoryorder`
            to *total ascending* or *total descending* if order
            should be determined by the numerical order of the
            values. Similarly, the order can be determined by the
            min, max, sum, mean, geometric mean or median of all
            the values.
        color
            Sets default for all colors associated with this axis
            all at once: line, font, tick, and grid colors. Grid
            color is lightened by blending this with the plot
            background Individual pieces can override this.
        constrain
            If this axis needs to be compressed (either due to its
            own `scaleanchor` and `scaleratio` or those of the
            other axis), determines how that happens: by increasing
            the "range", or by decreasing the "domain". Default is
            "domain" for axes containing image traces, "range"
            otherwise.
        constraintoward
            If this axis needs to be compressed (either due to its
            own `scaleanchor` and `scaleratio` or those of the
            other axis), determines which direction we push the
            originally specified plot area. Options are "left",
            "center" (default), and "right" for x axes, and "top",
            "middle" (default), and "bottom" for y axes.
        dividercolor
            Sets the color of the dividers Only has an effect on
            "multicategory" axes.
        dividerwidth
            Sets the width (in px) of the dividers Only has an
            effect on "multicategory" axes.
        domain
            Sets the domain of this axis (in plot fraction).
        dtick
            Sets the step in-between ticks on this axis. Use with
            `tick0`. Must be a positive number, or special strings
            available to "log" and "date" axes. If the axis `type`
            is "log", then ticks are set every 10^(n*dtick) where n
            is the tick number. For example, to set a tick mark at
            1, 10, 100, 1000, ... set dtick to 1. To set tick marks
            at 1, 100, 10000, ... set dtick to 2. To set tick marks
            at 1, 5, 25, 125, 625, 3125, ... set dtick to
            log_10(5), or 0.69897000433. "log" has several special
            values; "L<f>", where `f` is a positive number, gives
            ticks linearly spaced in value (but not position). For
            example `tick0` = 0.1, `dtick` = "L0.5" will put ticks
            at 0.1, 0.6, 1.1, 1.6 etc. To show powers of 10 plus
            small digits between, use "D1" (all digits) or "D2"
            (only 2 and 5). `tick0` is ignored for "D1" and "D2".
            If the axis `type` is "date", then you must convert the
            time to milliseconds. For example, to set the interval
            between ticks to one day, set `dtick` to 86400000.0.
            "date" also has special values "M<n>" gives ticks
            spaced by a number of months. `n` must be a positive
            integer. To set ticks on the 15th of every third month,
            set `tick0` to "2000-01-15" and `dtick` to "M3". To set
            ticks every 4 years, set `dtick` to "M48"
        exponentformat
            Determines a formatting rule for the tick exponents.
            For example, consider the number 1,000,000,000. If
            "none", it appears as 1,000,000,000. If "e", 1e+9. If
            "E", 1E+9. If "power", 1x10^9 (with 9 in a super
            script). If "SI", 1G. If "B", 1B.
        fixedrange
            Determines whether or not this axis is zoom-able. If
            true, then zoom is disabled.
        gridcolor
            Sets the color of the grid lines.
        griddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        gridwidth
            Sets the width (in px) of the grid lines.
        hoverformat
            Sets the hover text formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display "09~15~23.46"
        insiderange
            Could be used to set the desired inside range of this
            axis (excluding the labels) when `ticklabelposition` of
            the anchored axis has "inside". Not implemented for
            axes with `type` "log". This would be ignored when
            `range` is provided.
        labelalias
            Replacement text for specific tick or hover labels. For
            example using {US: 'USA', CA: 'Canada'} changes US to
            USA and CA to Canada. The labels we would have shown
            must match the keys exactly, after adding any
            tickprefix or ticksuffix. For negative numbers the
            minus sign symbol used (U+2212) is wider than the
            regular ascii dash. That means you need to use −1
            instead of -1. labelalias can be used with any axis
            type, and both keys (if needed) and values (if desired)
            can include html-like tags or MathJax.
        layer
            Sets the layer on which this axis is displayed. If
            *above traces*, this axis is displayed above all the
            subplot's traces If *below traces*, this axis is
            displayed below all the subplot's traces, but above the
            grid lines. Useful when used together with scatter-like
            traces with `cliponaxis` set to False to show markers
            and/or text nodes above this axis.
        linecolor
            Sets the axis line color.
        linewidth
            Sets the width (in px) of the axis line.
        matches
            If set to another axis id (e.g. `x2`, `y`), the range
            of this axis will match the range of the corresponding
            axis in data-coordinates space. Moreover, matching axes
            share auto-range values, category lists and histogram
            auto-bins. Note that setting axes simultaneously in
            both a `scaleanchor` and a `matches` constraint is
            currently forbidden. Moreover, note that matching axes
            must have the same `type`.
        maxallowed
            Determines the maximum range of this axis.
        minallowed
            Determines the minimum range of this axis.
        minexponent
            Hide SI prefix for 10^n if |n| is below this number.
            This only has an effect when `tickformat` is "SI" or
            "B".
        minor
            :class:`plotly.graph_objects.layout.xaxis.Minor`
            instance or dict with compatible properties
        mirror
            Determines if the axis lines or/and ticks are mirrored
            to the opposite side of the plotting area. If True, the
            axis lines are mirrored. If "ticks", the axis lines and
            ticks are mirrored. If False, mirroring is disable. If
            "all", axis lines are mirrored on all shared-axes
            subplots. If "allticks", axis lines and ticks are
            mirrored on all shared-axes subplots.
        nticks
            Specifies the maximum number of ticks for the
            particular axis. The actual number of ticks will be
            chosen automatically to be less than or equal to
            `nticks`. Has an effect only if `tickmode` is set to
            "auto".
        overlaying
            If set a same-letter axis id, this axis is overlaid on
            top of the corresponding same-letter axis, with traces
            and axes visible for both axes. If False, this axis
            does not overlay any same-letter axes. In this case,
            for axes with overlapping domains only the highest-
            numbered axis will be visible.
        position
            Sets the position of this axis in the plotting space
            (in normalized coordinates). Only has an effect if
            `anchor` is set to "free".
        range
            Sets the range of this axis. If the axis `type` is
            "log", then you must take the log of your desired range
            (e.g. to set the range from 1 to 100, set the range
            from 0 to 2). If the axis `type` is "date", it should
            be date strings, like date data, though Date objects
            and unix milliseconds will be accepted and converted to
            strings. If the axis `type` is "category", it should be
            numbers, using the scale where each category is
            assigned a serial number from zero in the order it
            appears. Leaving either or both elements `null` impacts
            the default `autorange`.
        rangebreaks
            A tuple of
            :class:`plotly.graph_objects.layout.xaxis.Rangebreak`
            instances or dicts with compatible properties
        rangebreakdefaults
            When used in a template (as
            layout.template.layout.xaxis.rangebreakdefaults), sets
            the default property values to use for elements of
            layout.xaxis.rangebreaks
        rangemode
            If "normal", the range is computed in relation to the
            extrema of the input data. If *tozero*`, the range
            extends to 0, regardless of the input data If
            "nonnegative", the range is non-negative, regardless of
            the input data. Applies only to linear axes.
        rangeselector
            :class:`plotly.graph_objects.layout.xaxis.Rangeselector
            ` instance or dict with compatible properties
        rangeslider
            :class:`plotly.graph_objects.layout.xaxis.Rangeslider`
            instance or dict with compatible properties
        scaleanchor
            If set to another axis id (e.g. `x2`, `y`), the range
            of this axis changes together with the range of the
            corresponding axis such that the scale of pixels per
            unit is in a constant ratio. Both axes are still
            zoomable, but when you zoom one, the other will zoom
            the same amount, keeping a fixed midpoint. `constrain`
            and `constraintoward` determine how we enforce the
            constraint. You can chain these, ie `yaxis:
            {scaleanchor: *x*}, xaxis2: {scaleanchor: *y*}` but you
            can only link axes of the same `type`. The linked axis
            can have the opposite letter (to constrain the aspect
            ratio) or the same letter (to match scales across
            subplots). Loops (`yaxis: {scaleanchor: *x*}, xaxis:
            {scaleanchor: *y*}` or longer) are redundant and the
            last constraint encountered will be ignored to avoid
            possible inconsistent constraints via `scaleratio`.
            Note that setting axes simultaneously in both a
            `scaleanchor` and a `matches` constraint is currently
            forbidden. Setting `false` allows to remove a default
            constraint (occasionally, you may need to prevent a
            default `scaleanchor` constraint from being applied,
            eg. when having an image trace `yaxis: {scaleanchor:
            "x"}` is set automatically in order for pixels to be
            rendered as squares, setting `yaxis: {scaleanchor:
            false}` allows to remove the constraint).
        scaleratio
            If this axis is linked to another by `scaleanchor`,
            this determines the pixel to unit scale ratio. For
            example, if this value is 10, then every unit on this
            axis spans 10 times the number of pixels as a unit on
            the linked axis. Use this for example to create an
            elevation profile where the vertical scale is
            exaggerated a fixed amount with respect to the
            horizontal.
        separatethousands
            If "true", even 4-digit integers are separated
        showdividers
            Determines whether or not a dividers are drawn between
            the category levels of this axis. Only has an effect on
            "multicategory" axes.
        showexponent
            If "all", all exponents are shown besides their
            significands. If "first", only the exponent of the
            first tick is shown. If "last", only the exponent of
            the last tick is shown. If "none", no exponents appear.
        showgrid
            Determines whether or not grid lines are drawn. If
            True, the grid lines are drawn at every tick mark.
        showline
            Determines whether or not a line bounding this axis is
            drawn.
        showspikes
            Determines whether or not spikes (aka droplines) are
            drawn for this axis. Note: This only takes affect when
            hovermode = closest
        showticklabels
            Determines whether or not the tick labels are drawn.
        showtickprefix
            If "all", all tick labels are displayed with a prefix.
            If "first", only the first tick is displayed with a
            prefix. If "last", only the last tick is displayed with
            a suffix. If "none", tick prefixes are hidden.
        showticksuffix
            Same as `showtickprefix` but for tick suffixes.
        side
            Determines whether a x (y) axis is positioned at the
            "bottom" ("left") or "top" ("right") of the plotting
            area.
        spikecolor
            Sets the spike color. If undefined, will use the series
            color
        spikedash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        spikemode
            Determines the drawing mode for the spike line If
            "toaxis", the line is drawn from the data point to the
            axis the  series is plotted on. If "across", the line
            is drawn across the entire plot area, and supercedes
            "toaxis". If "marker", then a marker dot is drawn on
            the axis the series is plotted on
        spikesnap
            Determines whether spikelines are stuck to the cursor
            or to the closest datapoints.
        spikethickness
            Sets the width (in px) of the zero line.
        tick0
            Sets the placement of the first tick on this axis. Use
            with `dtick`. If the axis `type` is "log", then you
            must take the log of your starting tick (e.g. to set
            the starting tick to 100, set the `tick0` to 2) except
            when `dtick`=*L<f>* (see `dtick` for more info). If the
            axis `type` is "date", it should be a date string, like
            date data. If the axis `type` is "category", it should
            be a number, using the scale where each category is
            assigned a serial number from zero in the order it
            appears.
        tickangle
            Sets the angle of the tick labels with respect to the
            horizontal. For example, a `tickangle` of -90 draws the
            tick labels vertically.
        tickcolor
            Sets the tick color.
        tickfont
            Sets the tick font.
        tickformat
            Sets the tick label formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display "09~15~23.46"
        tickformatstops
            A tuple of :class:`plotly.graph_objects.layout.xaxis.Ti
            ckformatstop` instances or dicts with compatible
            properties
        tickformatstopdefaults
            When used in a template (as
            layout.template.layout.xaxis.tickformatstopdefaults),
            sets the default property values to use for elements of
            layout.xaxis.tickformatstops
        ticklabelindex
            Only for axes with `type` "date" or "linear". Instead
            of drawing the major tick label, draw the label for the
            minor tick that is n positions away from the major
            tick. E.g. to always draw the label for the minor tick
            before each major tick, choose `ticklabelindex` -1.
            This is useful for date axes with `ticklabelmode`
            "period" if you want to label the period that ends with
            each major tick instead of the period that begins
            there.
        ticklabelindexsrc
            Sets the source reference on Chart Studio Cloud for
            `ticklabelindex`.
        ticklabelmode
            Determines where tick labels are drawn with respect to
            their corresponding ticks and grid lines. Only has an
            effect for axes of `type` "date" When set to "period",
            tick labels are drawn in the middle of the period
            between ticks.
        ticklabeloverflow
            Determines how we handle tick labels that would
            overflow either the graph div or the domain of the
            axis. The default value for inside tick labels is *hide
            past domain*. Otherwise on "category" and
            "multicategory" axes the default is "allow". In other
            cases the default is *hide past div*.
        ticklabelposition
            Determines where tick labels are drawn with respect to
            the axis Please note that top or bottom has no effect
            on x axes or when `ticklabelmode` is set to "period".
            Similarly left or right has no effect on y axes or when
            `ticklabelmode` is set to "period". Has no effect on
            "multicategory" axes or when `tickson` is set to
            "boundaries". When used on axes linked by `matches` or
            `scaleanchor`, no extra padding for inside labels would
            be added by autorange, so that the scales could match.
        ticklabelshift
            Shifts the tick labels by the specified number of
            pixels in parallel to the axis. Positive values move
            the labels in the positive direction of the axis.
        ticklabelstandoff
            Sets the standoff distance (in px) between the axis
            tick labels and their default position. A positive
            `ticklabelstandoff` moves the labels farther away from
            the plot area if `ticklabelposition` is "outside", and
            deeper into the plot area if `ticklabelposition` is
            "inside". A negative `ticklabelstandoff` works in the
            opposite direction, moving outside ticks towards the
            plot area and inside ticks towards the outside. If the
            negative value is large enough, inside ticks can even
            end up outside and vice versa.
        ticklabelstep
            Sets the spacing between tick labels as compared to the
            spacing between ticks. A value of 1 (default) means
            each tick gets a label. A value of 2 means shows every
            2nd label. A larger value n means only every nth tick
            is labeled. `tick0` determines which labels are shown.
            Not implemented for axes with `type` "log" or
            "multicategory", or when `tickmode` is "array".
        ticklen
            Sets the tick length (in px).
        tickmode
            Sets the tick mode for this axis. If "auto", the number
            of ticks is set via `nticks`. If "linear", the
            placement of the ticks is determined by a starting
            position `tick0` and a tick step `dtick` ("linear" is
            the default value if `tick0` and `dtick` are provided).
            If "array", the placement of the ticks is set via
            `tickvals` and the tick text is `ticktext`. ("array" is
            the default value if `tickvals` is provided). If
            "sync", the number of ticks will sync with the
            overlayed axis set by `overlaying` property.
        tickprefix
            Sets a tick label prefix.
        ticks
            Determines whether ticks are drawn or not. If "", this
            axis' ticks are not drawn. If "outside" ("inside"),
            this axis' are drawn outside (inside) the axis lines.
        tickson
            Determines where ticks and grid lines are drawn with
            respect to their corresponding tick labels. Only has an
            effect for axes of `type` "category" or
            "multicategory". When set to "boundaries", ticks and
            grid lines are drawn half a category to the left/bottom
            of labels.
        ticksuffix
            Sets a tick label suffix.
        ticktext
            Sets the text displayed at the ticks position via
            `tickvals`. Only has an effect if `tickmode` is set to
            "array". Used with `tickvals`.
        ticktextsrc
            Sets the source reference on Chart Studio Cloud for
            `ticktext`.
        tickvals
            Sets the values at which ticks on this axis appear.
            Only has an effect if `tickmode` is set to "array".
            Used with `ticktext`.
        tickvalssrc
            Sets the source reference on Chart Studio Cloud for
            `tickvals`.
        tickwidth
            Sets the tick width (in px).
        title
            :class:`plotly.graph_objects.layout.xaxis.Title`
            instance or dict with compatible properties
        titlefont
            Deprecated: Please use layout.xaxis.title.font instead.
            Sets this axis' title font. Note that the title's font
            used to be customized by the now deprecated `titlefont`
            attribute.
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.
        uirevision
            Controls persistence of user-driven changes in axis
            `range`, `autorange`, and `title` if in `editable:
            true` configuration. Defaults to `layout.uirevision`.
        visible
            A single toggle to hide the axis while preserving
            interaction like dragging. Default is true when a
            cheater plot is present on the axis, otherwise false
        zeroline
            Determines whether or not a line is drawn at along the
            0 value of this axis. If True, the zero line is drawn
            on top of the grid lines.
        zerolinecolor
            Sets the line color of the zero line.
        zerolinewidth
            Sets the width (in px) of the zero line.

        Returns
        -------
        XAxis
        """
        super(XAxis, self).__init__("xaxis")

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
The first argument to the plotly.graph_objs.layout.XAxis
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.XAxis`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("anchor", None)
        _v = anchor if anchor is not None else _v
        if _v is not None:
            self["anchor"] = _v
        _v = arg.pop("automargin", None)
        _v = automargin if automargin is not None else _v
        if _v is not None:
            self["automargin"] = _v
        _v = arg.pop("autorange", None)
        _v = autorange if autorange is not None else _v
        if _v is not None:
            self["autorange"] = _v
        _v = arg.pop("autorangeoptions", None)
        _v = autorangeoptions if autorangeoptions is not None else _v
        if _v is not None:
            self["autorangeoptions"] = _v
        _v = arg.pop("autotickangles", None)
        _v = autotickangles if autotickangles is not None else _v
        if _v is not None:
            self["autotickangles"] = _v
        _v = arg.pop("autotypenumbers", None)
        _v = autotypenumbers if autotypenumbers is not None else _v
        if _v is not None:
            self["autotypenumbers"] = _v
        _v = arg.pop("calendar", None)
        _v = calendar if calendar is not None else _v
        if _v is not None:
            self["calendar"] = _v
        _v = arg.pop("categoryarray", None)
        _v = categoryarray if categoryarray is not None else _v
        if _v is not None:
            self["categoryarray"] = _v
        _v = arg.pop("categoryarraysrc", None)
        _v = categoryarraysrc if categoryarraysrc is not None else _v
        if _v is not None:
            self["categoryarraysrc"] = _v
        _v = arg.pop("categoryorder", None)
        _v = categoryorder if categoryorder is not None else _v
        if _v is not None:
            self["categoryorder"] = _v
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("constrain", None)
        _v = constrain if constrain is not None else _v
        if _v is not None:
            self["constrain"] = _v
        _v = arg.pop("constraintoward", None)
        _v = constraintoward if constraintoward is not None else _v
        if _v is not None:
            self["constraintoward"] = _v
        _v = arg.pop("dividercolor", None)
        _v = dividercolor if dividercolor is not None else _v
        if _v is not None:
            self["dividercolor"] = _v
        _v = arg.pop("dividerwidth", None)
        _v = dividerwidth if dividerwidth is not None else _v
        if _v is not None:
            self["dividerwidth"] = _v
        _v = arg.pop("domain", None)
        _v = domain if domain is not None else _v
        if _v is not None:
            self["domain"] = _v
        _v = arg.pop("dtick", None)
        _v = dtick if dtick is not None else _v
        if _v is not None:
            self["dtick"] = _v
        _v = arg.pop("exponentformat", None)
        _v = exponentformat if exponentformat is not None else _v
        if _v is not None:
            self["exponentformat"] = _v
        _v = arg.pop("fixedrange", None)
        _v = fixedrange if fixedrange is not None else _v
        if _v is not None:
            self["fixedrange"] = _v
        _v = arg.pop("gridcolor", None)
        _v = gridcolor if gridcolor is not None else _v
        if _v is not None:
            self["gridcolor"] = _v
        _v = arg.pop("griddash", None)
        _v = griddash if griddash is not None else _v
        if _v is not None:
            self["griddash"] = _v
        _v = arg.pop("gridwidth", None)
        _v = gridwidth if gridwidth is not None else _v
        if _v is not None:
            self["gridwidth"] = _v
        _v = arg.pop("hoverformat", None)
        _v = hoverformat if hoverformat is not None else _v
        if _v is not None:
            self["hoverformat"] = _v
        _v = arg.pop("insiderange", None)
        _v = insiderange if insiderange is not None else _v
        if _v is not None:
            self["insiderange"] = _v
        _v = arg.pop("labelalias", None)
        _v = labelalias if labelalias is not None else _v
        if _v is not None:
            self["labelalias"] = _v
        _v = arg.pop("layer", None)
        _v = layer if layer is not None else _v
        if _v is not None:
            self["layer"] = _v
        _v = arg.pop("linecolor", None)
        _v = linecolor if linecolor is not None else _v
        if _v is not None:
            self["linecolor"] = _v
        _v = arg.pop("linewidth", None)
        _v = linewidth if linewidth is not None else _v
        if _v is not None:
            self["linewidth"] = _v
        _v = arg.pop("matches", None)
        _v = matches if matches is not None else _v
        if _v is not None:
            self["matches"] = _v
        _v = arg.pop("maxallowed", None)
        _v = maxallowed if maxallowed is not None else _v
        if _v is not None:
            self["maxallowed"] = _v
        _v = arg.pop("minallowed", None)
        _v = minallowed if minallowed is not None else _v
        if _v is not None:
            self["minallowed"] = _v
        _v = arg.pop("minexponent", None)
        _v = minexponent if minexponent is not None else _v
        if _v is not None:
            self["minexponent"] = _v
        _v = arg.pop("minor", None)
        _v = minor if minor is not None else _v
        if _v is not None:
            self["minor"] = _v
        _v = arg.pop("mirror", None)
        _v = mirror if mirror is not None else _v
        if _v is not None:
            self["mirror"] = _v
        _v = arg.pop("nticks", None)
        _v = nticks if nticks is not None else _v
        if _v is not None:
            self["nticks"] = _v
        _v = arg.pop("overlaying", None)
        _v = overlaying if overlaying is not None else _v
        if _v is not None:
            self["overlaying"] = _v
        _v = arg.pop("position", None)
        _v = position if position is not None else _v
        if _v is not None:
            self["position"] = _v
        _v = arg.pop("range", None)
        _v = range if range is not None else _v
        if _v is not None:
            self["range"] = _v
        _v = arg.pop("rangebreaks", None)
        _v = rangebreaks if rangebreaks is not None else _v
        if _v is not None:
            self["rangebreaks"] = _v
        _v = arg.pop("rangebreakdefaults", None)
        _v = rangebreakdefaults if rangebreakdefaults is not None else _v
        if _v is not None:
            self["rangebreakdefaults"] = _v
        _v = arg.pop("rangemode", None)
        _v = rangemode if rangemode is not None else _v
        if _v is not None:
            self["rangemode"] = _v
        _v = arg.pop("rangeselector", None)
        _v = rangeselector if rangeselector is not None else _v
        if _v is not None:
            self["rangeselector"] = _v
        _v = arg.pop("rangeslider", None)
        _v = rangeslider if rangeslider is not None else _v
        if _v is not None:
            self["rangeslider"] = _v
        _v = arg.pop("scaleanchor", None)
        _v = scaleanchor if scaleanchor is not None else _v
        if _v is not None:
            self["scaleanchor"] = _v
        _v = arg.pop("scaleratio", None)
        _v = scaleratio if scaleratio is not None else _v
        if _v is not None:
            self["scaleratio"] = _v
        _v = arg.pop("separatethousands", None)
        _v = separatethousands if separatethousands is not None else _v
        if _v is not None:
            self["separatethousands"] = _v
        _v = arg.pop("showdividers", None)
        _v = showdividers if showdividers is not None else _v
        if _v is not None:
            self["showdividers"] = _v
        _v = arg.pop("showexponent", None)
        _v = showexponent if showexponent is not None else _v
        if _v is not None:
            self["showexponent"] = _v
        _v = arg.pop("showgrid", None)
        _v = showgrid if showgrid is not None else _v
        if _v is not None:
            self["showgrid"] = _v
        _v = arg.pop("showline", None)
        _v = showline if showline is not None else _v
        if _v is not None:
            self["showline"] = _v
        _v = arg.pop("showspikes", None)
        _v = showspikes if showspikes is not None else _v
        if _v is not None:
            self["showspikes"] = _v
        _v = arg.pop("showticklabels", None)
        _v = showticklabels if showticklabels is not None else _v
        if _v is not None:
            self["showticklabels"] = _v
        _v = arg.pop("showtickprefix", None)
        _v = showtickprefix if showtickprefix is not None else _v
        if _v is not None:
            self["showtickprefix"] = _v
        _v = arg.pop("showticksuffix", None)
        _v = showticksuffix if showticksuffix is not None else _v
        if _v is not None:
            self["showticksuffix"] = _v
        _v = arg.pop("side", None)
        _v = side if side is not None else _v
        if _v is not None:
            self["side"] = _v
        _v = arg.pop("spikecolor", None)
        _v = spikecolor if spikecolor is not None else _v
        if _v is not None:
            self["spikecolor"] = _v
        _v = arg.pop("spikedash", None)
        _v = spikedash if spikedash is not None else _v
        if _v is not None:
            self["spikedash"] = _v
        _v = arg.pop("spikemode", None)
        _v = spikemode if spikemode is not None else _v
        if _v is not None:
            self["spikemode"] = _v
        _v = arg.pop("spikesnap", None)
        _v = spikesnap if spikesnap is not None else _v
        if _v is not None:
            self["spikesnap"] = _v
        _v = arg.pop("spikethickness", None)
        _v = spikethickness if spikethickness is not None else _v
        if _v is not None:
            self["spikethickness"] = _v
        _v = arg.pop("tick0", None)
        _v = tick0 if tick0 is not None else _v
        if _v is not None:
            self["tick0"] = _v
        _v = arg.pop("tickangle", None)
        _v = tickangle if tickangle is not None else _v
        if _v is not None:
            self["tickangle"] = _v
        _v = arg.pop("tickcolor", None)
        _v = tickcolor if tickcolor is not None else _v
        if _v is not None:
            self["tickcolor"] = _v
        _v = arg.pop("tickfont", None)
        _v = tickfont if tickfont is not None else _v
        if _v is not None:
            self["tickfont"] = _v
        _v = arg.pop("tickformat", None)
        _v = tickformat if tickformat is not None else _v
        if _v is not None:
            self["tickformat"] = _v
        _v = arg.pop("tickformatstops", None)
        _v = tickformatstops if tickformatstops is not None else _v
        if _v is not None:
            self["tickformatstops"] = _v
        _v = arg.pop("tickformatstopdefaults", None)
        _v = tickformatstopdefaults if tickformatstopdefaults is not None else _v
        if _v is not None:
            self["tickformatstopdefaults"] = _v
        _v = arg.pop("ticklabelindex", None)
        _v = ticklabelindex if ticklabelindex is not None else _v
        if _v is not None:
            self["ticklabelindex"] = _v
        _v = arg.pop("ticklabelindexsrc", None)
        _v = ticklabelindexsrc if ticklabelindexsrc is not None else _v
        if _v is not None:
            self["ticklabelindexsrc"] = _v
        _v = arg.pop("ticklabelmode", None)
        _v = ticklabelmode if ticklabelmode is not None else _v
        if _v is not None:
            self["ticklabelmode"] = _v
        _v = arg.pop("ticklabeloverflow", None)
        _v = ticklabeloverflow if ticklabeloverflow is not None else _v
        if _v is not None:
            self["ticklabeloverflow"] = _v
        _v = arg.pop("ticklabelposition", None)
        _v = ticklabelposition if ticklabelposition is not None else _v
        if _v is not None:
            self["ticklabelposition"] = _v
        _v = arg.pop("ticklabelshift", None)
        _v = ticklabelshift if ticklabelshift is not None else _v
        if _v is not None:
            self["ticklabelshift"] = _v
        _v = arg.pop("ticklabelstandoff", None)
        _v = ticklabelstandoff if ticklabelstandoff is not None else _v
        if _v is not None:
            self["ticklabelstandoff"] = _v
        _v = arg.pop("ticklabelstep", None)
        _v = ticklabelstep if ticklabelstep is not None else _v
        if _v is not None:
            self["ticklabelstep"] = _v
        _v = arg.pop("ticklen", None)
        _v = ticklen if ticklen is not None else _v
        if _v is not None:
            self["ticklen"] = _v
        _v = arg.pop("tickmode", None)
        _v = tickmode if tickmode is not None else _v
        if _v is not None:
            self["tickmode"] = _v
        _v = arg.pop("tickprefix", None)
        _v = tickprefix if tickprefix is not None else _v
        if _v is not None:
            self["tickprefix"] = _v
        _v = arg.pop("ticks", None)
        _v = ticks if ticks is not None else _v
        if _v is not None:
            self["ticks"] = _v
        _v = arg.pop("tickson", None)
        _v = tickson if tickson is not None else _v
        if _v is not None:
            self["tickson"] = _v
        _v = arg.pop("ticksuffix", None)
        _v = ticksuffix if ticksuffix is not None else _v
        if _v is not None:
            self["ticksuffix"] = _v
        _v = arg.pop("ticktext", None)
        _v = ticktext if ticktext is not None else _v
        if _v is not None:
            self["ticktext"] = _v
        _v = arg.pop("ticktextsrc", None)
        _v = ticktextsrc if ticktextsrc is not None else _v
        if _v is not None:
            self["ticktextsrc"] = _v
        _v = arg.pop("tickvals", None)
        _v = tickvals if tickvals is not None else _v
        if _v is not None:
            self["tickvals"] = _v
        _v = arg.pop("tickvalssrc", None)
        _v = tickvalssrc if tickvalssrc is not None else _v
        if _v is not None:
            self["tickvalssrc"] = _v
        _v = arg.pop("tickwidth", None)
        _v = tickwidth if tickwidth is not None else _v
        if _v is not None:
            self["tickwidth"] = _v
        _v = arg.pop("title", None)
        _v = title if title is not None else _v
        if _v is not None:
            self["title"] = _v
        _v = arg.pop("titlefont", None)
        _v = titlefont if titlefont is not None else _v
        if _v is not None:
            self["titlefont"] = _v
        _v = arg.pop("type", None)
        _v = type if type is not None else _v
        if _v is not None:
            self["type"] = _v
        _v = arg.pop("uirevision", None)
        _v = uirevision if uirevision is not None else _v
        if _v is not None:
            self["uirevision"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v
        _v = arg.pop("zeroline", None)
        _v = zeroline if zeroline is not None else _v
        if _v is not None:
            self["zeroline"] = _v
        _v = arg.pop("zerolinecolor", None)
        _v = zerolinecolor if zerolinecolor is not None else _v
        if _v is not None:
            self["zerolinecolor"] = _v
        _v = arg.pop("zerolinewidth", None)
        _v = zerolinewidth if zerolinewidth is not None else _v
        if _v is not None:
            self["zerolinewidth"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
