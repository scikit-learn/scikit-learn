#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class XAxis(_BaseLayoutHierarchyType):
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
        "minorloglabels",
        "mirror",
        "modebardisable",
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
        "type",
        "uirevision",
        "unifiedhovertitle",
        "visible",
        "zeroline",
        "zerolinecolor",
        "zerolinelayer",
        "zerolinewidth",
    }

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

    @property
    def autorangeoptions(self):
        """
        The 'autorangeoptions' property is an instance of Autorangeoptions
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Autorangeoptions`
          - A dict of string/value properties that will be passed
            to the Autorangeoptions constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Autorangeoptions
        """
        return self["autorangeoptions"]

    @autorangeoptions.setter
    def autorangeoptions(self, val):
        self["autorangeoptions"] = val

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
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

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
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["dividercolor"]

    @dividercolor.setter
    def dividercolor(self, val):
        self["dividercolor"] = val

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

    @property
    def gridcolor(self):
        """
        Sets the color of the grid lines.

        The 'gridcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["gridcolor"]

    @gridcolor.setter
    def gridcolor(self, val):
        self["gridcolor"] = val

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

    @property
    def linecolor(self):
        """
        Sets the axis line color.

        The 'linecolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["linecolor"]

    @linecolor.setter
    def linecolor(self, val):
        self["linecolor"] = val

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

    @property
    def minor(self):
        """
        The 'minor' property is an instance of Minor
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Minor`
          - A dict of string/value properties that will be passed
            to the Minor constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Minor
        """
        return self["minor"]

    @minor.setter
    def minor(self, val):
        self["minor"] = val

    @property
    def minorloglabels(self):
        """
        Determines how minor log labels are displayed. If *small
        digits*, small digits i.e. 2 or 5 are displayed. If "complete",
        complete digits are displayed. If "none", no labels are
        displayed.

        The 'minorloglabels' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['small digits', 'complete', 'none']

        Returns
        -------
        Any
        """
        return self["minorloglabels"]

    @minorloglabels.setter
    def minorloglabels(self, val):
        self["minorloglabels"] = val

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

    @property
    def modebardisable(self):
        """
        Disables certain modebar buttons for this axis. "autoscale"
        disables the autoscale buttons, "zoominout" disables the zoom-
        in and zoom-out buttons.

        The 'modebardisable' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['autoscale', 'zoominout'] joined with '+' characters
            (e.g. 'autoscale+zoominout')
            OR exactly one of ['none'] (e.g. 'none')

        Returns
        -------
        Any
        """
        return self["modebardisable"]

    @modebardisable.setter
    def modebardisable(self, val):
        self["modebardisable"] = val

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

    @property
    def rangebreaks(self):
        """
        The 'rangebreaks' property is a tuple of instances of
        Rangebreak that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.xaxis.Rangebreak
          - A list or tuple of dicts of string/value properties that
            will be passed to the Rangebreak constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.xaxis.Rangebreak]
        """
        return self["rangebreaks"]

    @rangebreaks.setter
    def rangebreaks(self, val):
        self["rangebreaks"] = val

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

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Rangebreak
        """
        return self["rangebreakdefaults"]

    @rangebreakdefaults.setter
    def rangebreakdefaults(self, val):
        self["rangebreakdefaults"] = val

    @property
    def rangemode(self):
        """
        If "normal", the range is computed in relation to the extrema
        of the input data. If "tozero", the range extends to 0,
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

    @property
    def rangeselector(self):
        """
        The 'rangeselector' property is an instance of Rangeselector
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Rangeselector`
          - A dict of string/value properties that will be passed
            to the Rangeselector constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Rangeselector
        """
        return self["rangeselector"]

    @rangeselector.setter
    def rangeselector(self, val):
        self["rangeselector"] = val

    @property
    def rangeslider(self):
        """
        The 'rangeslider' property is an instance of Rangeslider
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Rangeslider`
          - A dict of string/value properties that will be passed
            to the Rangeslider constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Rangeslider
        """
        return self["rangeslider"]

    @rangeslider.setter
    def rangeslider(self, val):
        self["rangeslider"] = val

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

    @property
    def spikecolor(self):
        """
        Sets the spike color. If undefined, will use the series color

        The 'spikecolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["spikecolor"]

    @spikecolor.setter
    def spikecolor(self, val):
        self["spikecolor"] = val

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

    @property
    def tickcolor(self):
        """
        Sets the tick color.

        The 'tickcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["tickcolor"]

    @tickcolor.setter
    def tickcolor(self, val):
        self["tickcolor"] = val

    @property
    def tickfont(self):
        """
        Sets the tick font.

        The 'tickfont' property is an instance of Tickfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Tickfont`
          - A dict of string/value properties that will be passed
            to the Tickfont constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Tickfont
        """
        return self["tickfont"]

    @tickfont.setter
    def tickfont(self, val):
        self["tickfont"] = val

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

    @property
    def tickformatstops(self):
        """
        The 'tickformatstops' property is a tuple of instances of
        Tickformatstop that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.xaxis.Tickformatstop
          - A list or tuple of dicts of string/value properties that
            will be passed to the Tickformatstop constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.xaxis.Tickformatstop]
        """
        return self["tickformatstops"]

    @tickformatstops.setter
    def tickformatstops(self, val):
        self["tickformatstops"] = val

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

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Tickformatstop
        """
        return self["tickformatstopdefaults"]

    @tickformatstopdefaults.setter
    def tickformatstopdefaults(self, val):
        self["tickformatstopdefaults"] = val

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

    @property
    def ticklabelposition(self):
        """
        Determines where tick labels are drawn with respect to the
        axis. Please note that top or bottom has no effect on x axes or
        when `ticklabelmode` is set to "period" or when `tickson` is
        set to "boundaries". Similarly, left or right has no effect on
        y axes or when `ticklabelmode` is set to "period" or when
        `tickson` is set to "boundaries". Has no effect on
        "multicategory" axes. When used on axes linked by `matches` or
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

    @property
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Title
        """
        return self["title"]

    @title.setter
    def title(self, val):
        self["title"] = val

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

    @property
    def unifiedhovertitle(self):
        """
        The 'unifiedhovertitle' property is an instance of Unifiedhovertitle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.Unifiedhovertitle`
          - A dict of string/value properties that will be passed
            to the Unifiedhovertitle constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.Unifiedhovertitle
        """
        return self["unifiedhovertitle"]

    @unifiedhovertitle.setter
    def unifiedhovertitle(self, val):
        self["unifiedhovertitle"] = val

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

    @property
    def zerolinecolor(self):
        """
        Sets the line color of the zero line.

        The 'zerolinecolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["zerolinecolor"]

    @zerolinecolor.setter
    def zerolinecolor(self, val):
        self["zerolinecolor"] = val

    @property
    def zerolinelayer(self):
        """
        Sets the layer on which this zeroline is displayed. If *above
        traces*, this zeroline is displayed above all the subplot's
        traces If *below traces*, this zeroline is displayed below all
        the subplot's traces, but above the grid lines. Limitation:
        "zerolinelayer" currently has no effect if the "zorder"
        property is set on any trace.

        The 'zerolinelayer' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['above traces', 'below traces']

        Returns
        -------
        Any
        """
        return self["zerolinelayer"]

    @zerolinelayer.setter
    def zerolinelayer(self, val):
        self["zerolinelayer"] = val

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
        minorloglabels
            Determines how minor log labels are displayed. If
            *small digits*, small digits i.e. 2 or 5 are displayed.
            If "complete", complete digits are displayed. If
            "none", no labels are displayed.
        mirror
            Determines if the axis lines or/and ticks are mirrored
            to the opposite side of the plotting area. If True, the
            axis lines are mirrored. If "ticks", the axis lines and
            ticks are mirrored. If False, mirroring is disable. If
            "all", axis lines are mirrored on all shared-axes
            subplots. If "allticks", axis lines and ticks are
            mirrored on all shared-axes subplots.
        modebardisable
            Disables certain modebar buttons for this axis.
            "autoscale" disables the autoscale buttons, "zoominout"
            disables the zoom-in and zoom-out buttons.
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
            extrema of the input data. If "tozero", the range
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
            the axis. Please note that top or bottom has no effect
            on x axes or when `ticklabelmode` is set to "period" or
            when `tickson` is set to "boundaries". Similarly, left
            or right has no effect on y axes or when
            `ticklabelmode` is set to "period" or when `tickson` is
            set to "boundaries". Has no effect on "multicategory"
            axes. When used on axes linked by `matches` or
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
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.
        uirevision
            Controls persistence of user-driven changes in axis
            `range`, `autorange`, and `title` if in `editable:
            true` configuration. Defaults to `layout.uirevision`.
        unifiedhovertitle
            :class:`plotly.graph_objects.layout.xaxis.Unifiedhovert
            itle` instance or dict with compatible properties
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
        zerolinelayer
            Sets the layer on which this zeroline is displayed. If
            *above traces*, this zeroline is displayed above all
            the subplot's traces If *below traces*, this zeroline
            is displayed below all the subplot's traces, but above
            the grid lines. Limitation: "zerolinelayer" currently
            has no effect if the "zorder" property is set on any
            trace.
        zerolinewidth
            Sets the width (in px) of the zero line.
        """

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
        minorloglabels=None,
        mirror=None,
        modebardisable=None,
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
        type=None,
        uirevision=None,
        unifiedhovertitle=None,
        visible=None,
        zeroline=None,
        zerolinecolor=None,
        zerolinelayer=None,
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
        minorloglabels
            Determines how minor log labels are displayed. If
            *small digits*, small digits i.e. 2 or 5 are displayed.
            If "complete", complete digits are displayed. If
            "none", no labels are displayed.
        mirror
            Determines if the axis lines or/and ticks are mirrored
            to the opposite side of the plotting area. If True, the
            axis lines are mirrored. If "ticks", the axis lines and
            ticks are mirrored. If False, mirroring is disable. If
            "all", axis lines are mirrored on all shared-axes
            subplots. If "allticks", axis lines and ticks are
            mirrored on all shared-axes subplots.
        modebardisable
            Disables certain modebar buttons for this axis.
            "autoscale" disables the autoscale buttons, "zoominout"
            disables the zoom-in and zoom-out buttons.
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
            extrema of the input data. If "tozero", the range
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
            the axis. Please note that top or bottom has no effect
            on x axes or when `ticklabelmode` is set to "period" or
            when `tickson` is set to "boundaries". Similarly, left
            or right has no effect on y axes or when
            `ticklabelmode` is set to "period" or when `tickson` is
            set to "boundaries". Has no effect on "multicategory"
            axes. When used on axes linked by `matches` or
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
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.
        uirevision
            Controls persistence of user-driven changes in axis
            `range`, `autorange`, and `title` if in `editable:
            true` configuration. Defaults to `layout.uirevision`.
        unifiedhovertitle
            :class:`plotly.graph_objects.layout.xaxis.Unifiedhovert
            itle` instance or dict with compatible properties
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
        zerolinelayer
            Sets the layer on which this zeroline is displayed. If
            *above traces*, this zeroline is displayed above all
            the subplot's traces If *below traces*, this zeroline
            is displayed below all the subplot's traces, but above
            the grid lines. Limitation: "zerolinelayer" currently
            has no effect if the "zorder" property is set on any
            trace.
        zerolinewidth
            Sets the width (in px) of the zero line.

        Returns
        -------
        XAxis
        """
        super().__init__("xaxis")
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
The first argument to the plotly.graph_objs.layout.XAxis
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.XAxis`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("anchor", arg, anchor)
        self._set_property("automargin", arg, automargin)
        self._set_property("autorange", arg, autorange)
        self._set_property("autorangeoptions", arg, autorangeoptions)
        self._set_property("autotickangles", arg, autotickangles)
        self._set_property("autotypenumbers", arg, autotypenumbers)
        self._set_property("calendar", arg, calendar)
        self._set_property("categoryarray", arg, categoryarray)
        self._set_property("categoryarraysrc", arg, categoryarraysrc)
        self._set_property("categoryorder", arg, categoryorder)
        self._set_property("color", arg, color)
        self._set_property("constrain", arg, constrain)
        self._set_property("constraintoward", arg, constraintoward)
        self._set_property("dividercolor", arg, dividercolor)
        self._set_property("dividerwidth", arg, dividerwidth)
        self._set_property("domain", arg, domain)
        self._set_property("dtick", arg, dtick)
        self._set_property("exponentformat", arg, exponentformat)
        self._set_property("fixedrange", arg, fixedrange)
        self._set_property("gridcolor", arg, gridcolor)
        self._set_property("griddash", arg, griddash)
        self._set_property("gridwidth", arg, gridwidth)
        self._set_property("hoverformat", arg, hoverformat)
        self._set_property("insiderange", arg, insiderange)
        self._set_property("labelalias", arg, labelalias)
        self._set_property("layer", arg, layer)
        self._set_property("linecolor", arg, linecolor)
        self._set_property("linewidth", arg, linewidth)
        self._set_property("matches", arg, matches)
        self._set_property("maxallowed", arg, maxallowed)
        self._set_property("minallowed", arg, minallowed)
        self._set_property("minexponent", arg, minexponent)
        self._set_property("minor", arg, minor)
        self._set_property("minorloglabels", arg, minorloglabels)
        self._set_property("mirror", arg, mirror)
        self._set_property("modebardisable", arg, modebardisable)
        self._set_property("nticks", arg, nticks)
        self._set_property("overlaying", arg, overlaying)
        self._set_property("position", arg, position)
        self._set_property("range", arg, range)
        self._set_property("rangebreaks", arg, rangebreaks)
        self._set_property("rangebreakdefaults", arg, rangebreakdefaults)
        self._set_property("rangemode", arg, rangemode)
        self._set_property("rangeselector", arg, rangeselector)
        self._set_property("rangeslider", arg, rangeslider)
        self._set_property("scaleanchor", arg, scaleanchor)
        self._set_property("scaleratio", arg, scaleratio)
        self._set_property("separatethousands", arg, separatethousands)
        self._set_property("showdividers", arg, showdividers)
        self._set_property("showexponent", arg, showexponent)
        self._set_property("showgrid", arg, showgrid)
        self._set_property("showline", arg, showline)
        self._set_property("showspikes", arg, showspikes)
        self._set_property("showticklabels", arg, showticklabels)
        self._set_property("showtickprefix", arg, showtickprefix)
        self._set_property("showticksuffix", arg, showticksuffix)
        self._set_property("side", arg, side)
        self._set_property("spikecolor", arg, spikecolor)
        self._set_property("spikedash", arg, spikedash)
        self._set_property("spikemode", arg, spikemode)
        self._set_property("spikesnap", arg, spikesnap)
        self._set_property("spikethickness", arg, spikethickness)
        self._set_property("tick0", arg, tick0)
        self._set_property("tickangle", arg, tickangle)
        self._set_property("tickcolor", arg, tickcolor)
        self._set_property("tickfont", arg, tickfont)
        self._set_property("tickformat", arg, tickformat)
        self._set_property("tickformatstops", arg, tickformatstops)
        self._set_property("tickformatstopdefaults", arg, tickformatstopdefaults)
        self._set_property("ticklabelindex", arg, ticklabelindex)
        self._set_property("ticklabelindexsrc", arg, ticklabelindexsrc)
        self._set_property("ticklabelmode", arg, ticklabelmode)
        self._set_property("ticklabeloverflow", arg, ticklabeloverflow)
        self._set_property("ticklabelposition", arg, ticklabelposition)
        self._set_property("ticklabelshift", arg, ticklabelshift)
        self._set_property("ticklabelstandoff", arg, ticklabelstandoff)
        self._set_property("ticklabelstep", arg, ticklabelstep)
        self._set_property("ticklen", arg, ticklen)
        self._set_property("tickmode", arg, tickmode)
        self._set_property("tickprefix", arg, tickprefix)
        self._set_property("ticks", arg, ticks)
        self._set_property("tickson", arg, tickson)
        self._set_property("ticksuffix", arg, ticksuffix)
        self._set_property("ticktext", arg, ticktext)
        self._set_property("ticktextsrc", arg, ticktextsrc)
        self._set_property("tickvals", arg, tickvals)
        self._set_property("tickvalssrc", arg, tickvalssrc)
        self._set_property("tickwidth", arg, tickwidth)
        self._set_property("title", arg, title)
        self._set_property("type", arg, type)
        self._set_property("uirevision", arg, uirevision)
        self._set_property("unifiedhovertitle", arg, unifiedhovertitle)
        self._set_property("visible", arg, visible)
        self._set_property("zeroline", arg, zeroline)
        self._set_property("zerolinecolor", arg, zerolinecolor)
        self._set_property("zerolinelayer", arg, zerolinelayer)
        self._set_property("zerolinewidth", arg, zerolinewidth)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
