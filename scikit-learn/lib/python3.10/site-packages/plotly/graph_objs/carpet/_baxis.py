#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Baxis(_BaseTraceHierarchyType):
    _parent_path_str = "carpet"
    _path_str = "carpet.baxis"
    _valid_props = {
        "arraydtick",
        "arraytick0",
        "autorange",
        "autotypenumbers",
        "categoryarray",
        "categoryarraysrc",
        "categoryorder",
        "cheatertype",
        "color",
        "dtick",
        "endline",
        "endlinecolor",
        "endlinewidth",
        "exponentformat",
        "fixedrange",
        "gridcolor",
        "griddash",
        "gridwidth",
        "labelalias",
        "labelpadding",
        "labelprefix",
        "labelsuffix",
        "linecolor",
        "linewidth",
        "minexponent",
        "minorgridcolor",
        "minorgridcount",
        "minorgriddash",
        "minorgridwidth",
        "nticks",
        "range",
        "rangemode",
        "separatethousands",
        "showexponent",
        "showgrid",
        "showline",
        "showticklabels",
        "showtickprefix",
        "showticksuffix",
        "smoothing",
        "startline",
        "startlinecolor",
        "startlinewidth",
        "tick0",
        "tickangle",
        "tickfont",
        "tickformat",
        "tickformatstopdefaults",
        "tickformatstops",
        "tickmode",
        "tickprefix",
        "ticksuffix",
        "ticktext",
        "ticktextsrc",
        "tickvals",
        "tickvalssrc",
        "title",
        "type",
    }

    @property
    def arraydtick(self):
        """
        The stride between grid lines along the axis

        The 'arraydtick' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["arraydtick"]

    @arraydtick.setter
    def arraydtick(self, val):
        self["arraydtick"] = val

    @property
    def arraytick0(self):
        """
        The starting index of grid lines along the axis

        The 'arraytick0' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["arraytick0"]

    @arraytick0.setter
    def arraytick0(self, val):
        self["arraytick0"] = val

    @property
    def autorange(self):
        """
        Determines whether or not the range of this axis is computed in
        relation to the input data. See `rangemode` for more info. If
        `range` is provided, then `autorange` is set to False.

        The 'autorange' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, False, 'reversed']

        Returns
        -------
        Any
        """
        return self["autorange"]

    @autorange.setter
    def autorange(self, val):
        self["autorange"] = val

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
        `categoryarray`.

        The 'categoryorder' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['trace', 'category ascending', 'category descending',
                'array']

        Returns
        -------
        Any
        """
        return self["categoryorder"]

    @categoryorder.setter
    def categoryorder(self, val):
        self["categoryorder"] = val

    @property
    def cheatertype(self):
        """
        The 'cheatertype' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['index', 'value']

        Returns
        -------
        Any
        """
        return self["cheatertype"]

    @cheatertype.setter
    def cheatertype(self, val):
        self["cheatertype"] = val

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
    def dtick(self):
        """
        The stride between grid lines along the axis

        The 'dtick' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["dtick"]

    @dtick.setter
    def dtick(self, val):
        self["dtick"] = val

    @property
    def endline(self):
        """
        Determines whether or not a line is drawn at along the final
        value of this axis. If True, the end line is drawn on top of
        the grid lines.

        The 'endline' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["endline"]

    @endline.setter
    def endline(self, val):
        self["endline"] = val

    @property
    def endlinecolor(self):
        """
        Sets the line color of the end line.

        The 'endlinecolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["endlinecolor"]

    @endlinecolor.setter
    def endlinecolor(self, val):
        self["endlinecolor"] = val

    @property
    def endlinewidth(self):
        """
        Sets the width (in px) of the end line.

        The 'endlinewidth' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["endlinewidth"]

    @endlinewidth.setter
    def endlinewidth(self, val):
        self["endlinewidth"] = val

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
        Sets the axis line color.

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
        Sets the width (in px) of the axis line.

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
    def labelpadding(self):
        """
        Extra padding between label and the axis

        The 'labelpadding' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)

        Returns
        -------
        int
        """
        return self["labelpadding"]

    @labelpadding.setter
    def labelpadding(self, val):
        self["labelpadding"] = val

    @property
    def labelprefix(self):
        """
        Sets a axis label prefix.

        The 'labelprefix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["labelprefix"]

    @labelprefix.setter
    def labelprefix(self, val):
        self["labelprefix"] = val

    @property
    def labelsuffix(self):
        """
        Sets a axis label suffix.

        The 'labelsuffix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["labelsuffix"]

    @labelsuffix.setter
    def labelsuffix(self, val):
        self["labelsuffix"] = val

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
    def minexponent(self):
        """
        Hide SI prefix for 10^n if |n| is below this number

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
    def minorgridcolor(self):
        """
        Sets the color of the grid lines.

        The 'minorgridcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["minorgridcolor"]

    @minorgridcolor.setter
    def minorgridcolor(self, val):
        self["minorgridcolor"] = val

    @property
    def minorgridcount(self):
        """
        Sets the number of minor grid ticks per major grid tick

        The 'minorgridcount' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["minorgridcount"]

    @minorgridcount.setter
    def minorgridcount(self, val):
        self["minorgridcount"] = val

    @property
    def minorgriddash(self):
        """
        Sets the dash style of lines. Set to a dash type string
        ("solid", "dot", "dash", "longdash", "dashdot", or
        "longdashdot") or a dash length list in px (eg
        "5px,10px,2px,2px").

        The 'minorgriddash' property is an enumeration that may be specified as:
          - One of the following dash styles:
                ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
          - A string containing a dash length list in pixels or percentages
                (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)

        Returns
        -------
        str
        """
        return self["minorgriddash"]

    @minorgriddash.setter
    def minorgriddash(self, val):
        self["minorgriddash"] = val

    @property
    def minorgridwidth(self):
        """
        Sets the width (in px) of the grid lines.

        The 'minorgridwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["minorgridwidth"]

    @minorgridwidth.setter
    def minorgridwidth(self, val):
        self["minorgridwidth"] = val

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
        If "normal", the range is computed in relation to the extrema
        of the input data. If "tozero", the range extends to 0,
        regardless of the input data If "nonnegative", the range is
        non-negative, regardless of the input data.

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
    def showticklabels(self):
        """
        Determines whether axis labels are drawn on the low side, the
        high side, both, or neither side of the axis.

        The 'showticklabels' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['start', 'end', 'both', 'none']

        Returns
        -------
        Any
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
    def smoothing(self):
        """
        The 'smoothing' property is a number and may be specified as:
          - An int or float in the interval [0, 1.3]

        Returns
        -------
        int|float
        """
        return self["smoothing"]

    @smoothing.setter
    def smoothing(self, val):
        self["smoothing"] = val

    @property
    def startline(self):
        """
        Determines whether or not a line is drawn at along the starting
        value of this axis. If True, the start line is drawn on top of
        the grid lines.

        The 'startline' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["startline"]

    @startline.setter
    def startline(self, val):
        self["startline"] = val

    @property
    def startlinecolor(self):
        """
        Sets the line color of the start line.

        The 'startlinecolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["startlinecolor"]

    @startlinecolor.setter
    def startlinecolor(self, val):
        self["startlinecolor"] = val

    @property
    def startlinewidth(self):
        """
        Sets the width (in px) of the start line.

        The 'startlinewidth' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["startlinewidth"]

    @startlinewidth.setter
    def startlinewidth(self, val):
        self["startlinewidth"] = val

    @property
    def tick0(self):
        """
        The starting index of grid lines along the axis

        The 'tick0' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
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
    def tickfont(self):
        """
        Sets the tick font.

        The 'tickfont' property is an instance of Tickfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.baxis.Tickfont`
          - A dict of string/value properties that will be passed
            to the Tickfont constructor

        Returns
        -------
        plotly.graph_objs.carpet.baxis.Tickfont
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
          - A list or tuple of instances of plotly.graph_objs.carpet.baxis.Tickformatstop
          - A list or tuple of dicts of string/value properties that
            will be passed to the Tickformatstop constructor

        Returns
        -------
        tuple[plotly.graph_objs.carpet.baxis.Tickformatstop]
        """
        return self["tickformatstops"]

    @tickformatstops.setter
    def tickformatstops(self, val):
        self["tickformatstops"] = val

    @property
    def tickformatstopdefaults(self):
        """
        When used in a template (as
        layout.template.data.carpet.baxis.tickformatstopdefaults), sets
        the default property values to use for elements of
        carpet.baxis.tickformatstops

        The 'tickformatstopdefaults' property is an instance of Tickformatstop
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.baxis.Tickformatstop`
          - A dict of string/value properties that will be passed
            to the Tickformatstop constructor

        Returns
        -------
        plotly.graph_objs.carpet.baxis.Tickformatstop
        """
        return self["tickformatstopdefaults"]

    @tickformatstopdefaults.setter
    def tickformatstopdefaults(self, val):
        self["tickformatstopdefaults"] = val

    @property
    def tickmode(self):
        """
        The 'tickmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['linear', 'array']

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
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.baxis.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

        Returns
        -------
        plotly.graph_objs.carpet.baxis.Title
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
                ['-', 'linear', 'date', 'category']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    @property
    def _prop_descriptions(self):
        return """\
        arraydtick
            The stride between grid lines along the axis
        arraytick0
            The starting index of grid lines along the axis
        autorange
            Determines whether or not the range of this axis is
            computed in relation to the input data. See `rangemode`
            for more info. If `range` is provided, then `autorange`
            is set to False.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. Defaults to
            layout.autotypenumbers.
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
            the categories in `categoryarray`.
        cheatertype

        color
            Sets default for all colors associated with this axis
            all at once: line, font, tick, and grid colors. Grid
            color is lightened by blending this with the plot
            background Individual pieces can override this.
        dtick
            The stride between grid lines along the axis
        endline
            Determines whether or not a line is drawn at along the
            final value of this axis. If True, the end line is
            drawn on top of the grid lines.
        endlinecolor
            Sets the line color of the end line.
        endlinewidth
            Sets the width (in px) of the end line.
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
            Sets the axis line color.
        griddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        gridwidth
            Sets the width (in px) of the axis line.
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
        labelpadding
            Extra padding between label and the axis
        labelprefix
            Sets a axis label prefix.
        labelsuffix
            Sets a axis label suffix.
        linecolor
            Sets the axis line color.
        linewidth
            Sets the width (in px) of the axis line.
        minexponent
            Hide SI prefix for 10^n if |n| is below this number
        minorgridcolor
            Sets the color of the grid lines.
        minorgridcount
            Sets the number of minor grid ticks per major grid tick
        minorgriddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        minorgridwidth
            Sets the width (in px) of the grid lines.
        nticks
            Specifies the maximum number of ticks for the
            particular axis. The actual number of ticks will be
            chosen automatically to be less than or equal to
            `nticks`. Has an effect only if `tickmode` is set to
            "auto".
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
            appears.
        rangemode
            If "normal", the range is computed in relation to the
            extrema of the input data. If "tozero", the range
            extends to 0, regardless of the input data If
            "nonnegative", the range is non-negative, regardless of
            the input data.
        separatethousands
            If "true", even 4-digit integers are separated
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
        showticklabels
            Determines whether axis labels are drawn on the low
            side, the high side, both, or neither side of the axis.
        showtickprefix
            If "all", all tick labels are displayed with a prefix.
            If "first", only the first tick is displayed with a
            prefix. If "last", only the last tick is displayed with
            a suffix. If "none", tick prefixes are hidden.
        showticksuffix
            Same as `showtickprefix` but for tick suffixes.
        smoothing

        startline
            Determines whether or not a line is drawn at along the
            starting value of this axis. If True, the start line is
            drawn on top of the grid lines.
        startlinecolor
            Sets the line color of the start line.
        startlinewidth
            Sets the width (in px) of the start line.
        tick0
            The starting index of grid lines along the axis
        tickangle
            Sets the angle of the tick labels with respect to the
            horizontal. For example, a `tickangle` of -90 draws the
            tick labels vertically.
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
            A tuple of :class:`plotly.graph_objects.carpet.baxis.Ti
            ckformatstop` instances or dicts with compatible
            properties
        tickformatstopdefaults
            When used in a template (as layout.template.data.carpet
            .baxis.tickformatstopdefaults), sets the default
            property values to use for elements of
            carpet.baxis.tickformatstops
        tickmode

        tickprefix
            Sets a tick label prefix.
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
        title
            :class:`plotly.graph_objects.carpet.baxis.Title`
            instance or dict with compatible properties
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.
        """

    def __init__(
        self,
        arg=None,
        arraydtick=None,
        arraytick0=None,
        autorange=None,
        autotypenumbers=None,
        categoryarray=None,
        categoryarraysrc=None,
        categoryorder=None,
        cheatertype=None,
        color=None,
        dtick=None,
        endline=None,
        endlinecolor=None,
        endlinewidth=None,
        exponentformat=None,
        fixedrange=None,
        gridcolor=None,
        griddash=None,
        gridwidth=None,
        labelalias=None,
        labelpadding=None,
        labelprefix=None,
        labelsuffix=None,
        linecolor=None,
        linewidth=None,
        minexponent=None,
        minorgridcolor=None,
        minorgridcount=None,
        minorgriddash=None,
        minorgridwidth=None,
        nticks=None,
        range=None,
        rangemode=None,
        separatethousands=None,
        showexponent=None,
        showgrid=None,
        showline=None,
        showticklabels=None,
        showtickprefix=None,
        showticksuffix=None,
        smoothing=None,
        startline=None,
        startlinecolor=None,
        startlinewidth=None,
        tick0=None,
        tickangle=None,
        tickfont=None,
        tickformat=None,
        tickformatstops=None,
        tickformatstopdefaults=None,
        tickmode=None,
        tickprefix=None,
        ticksuffix=None,
        ticktext=None,
        ticktextsrc=None,
        tickvals=None,
        tickvalssrc=None,
        title=None,
        type=None,
        **kwargs,
    ):
        """
        Construct a new Baxis object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.carpet.Baxis`
        arraydtick
            The stride between grid lines along the axis
        arraytick0
            The starting index of grid lines along the axis
        autorange
            Determines whether or not the range of this axis is
            computed in relation to the input data. See `rangemode`
            for more info. If `range` is provided, then `autorange`
            is set to False.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. Defaults to
            layout.autotypenumbers.
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
            the categories in `categoryarray`.
        cheatertype

        color
            Sets default for all colors associated with this axis
            all at once: line, font, tick, and grid colors. Grid
            color is lightened by blending this with the plot
            background Individual pieces can override this.
        dtick
            The stride between grid lines along the axis
        endline
            Determines whether or not a line is drawn at along the
            final value of this axis. If True, the end line is
            drawn on top of the grid lines.
        endlinecolor
            Sets the line color of the end line.
        endlinewidth
            Sets the width (in px) of the end line.
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
            Sets the axis line color.
        griddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        gridwidth
            Sets the width (in px) of the axis line.
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
        labelpadding
            Extra padding between label and the axis
        labelprefix
            Sets a axis label prefix.
        labelsuffix
            Sets a axis label suffix.
        linecolor
            Sets the axis line color.
        linewidth
            Sets the width (in px) of the axis line.
        minexponent
            Hide SI prefix for 10^n if |n| is below this number
        minorgridcolor
            Sets the color of the grid lines.
        minorgridcount
            Sets the number of minor grid ticks per major grid tick
        minorgriddash
            Sets the dash style of lines. Set to a dash type string
            ("solid", "dot", "dash", "longdash", "dashdot", or
            "longdashdot") or a dash length list in px (eg
            "5px,10px,2px,2px").
        minorgridwidth
            Sets the width (in px) of the grid lines.
        nticks
            Specifies the maximum number of ticks for the
            particular axis. The actual number of ticks will be
            chosen automatically to be less than or equal to
            `nticks`. Has an effect only if `tickmode` is set to
            "auto".
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
            appears.
        rangemode
            If "normal", the range is computed in relation to the
            extrema of the input data. If "tozero", the range
            extends to 0, regardless of the input data If
            "nonnegative", the range is non-negative, regardless of
            the input data.
        separatethousands
            If "true", even 4-digit integers are separated
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
        showticklabels
            Determines whether axis labels are drawn on the low
            side, the high side, both, or neither side of the axis.
        showtickprefix
            If "all", all tick labels are displayed with a prefix.
            If "first", only the first tick is displayed with a
            prefix. If "last", only the last tick is displayed with
            a suffix. If "none", tick prefixes are hidden.
        showticksuffix
            Same as `showtickprefix` but for tick suffixes.
        smoothing

        startline
            Determines whether or not a line is drawn at along the
            starting value of this axis. If True, the start line is
            drawn on top of the grid lines.
        startlinecolor
            Sets the line color of the start line.
        startlinewidth
            Sets the width (in px) of the start line.
        tick0
            The starting index of grid lines along the axis
        tickangle
            Sets the angle of the tick labels with respect to the
            horizontal. For example, a `tickangle` of -90 draws the
            tick labels vertically.
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
            A tuple of :class:`plotly.graph_objects.carpet.baxis.Ti
            ckformatstop` instances or dicts with compatible
            properties
        tickformatstopdefaults
            When used in a template (as layout.template.data.carpet
            .baxis.tickformatstopdefaults), sets the default
            property values to use for elements of
            carpet.baxis.tickformatstops
        tickmode

        tickprefix
            Sets a tick label prefix.
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
        title
            :class:`plotly.graph_objects.carpet.baxis.Title`
            instance or dict with compatible properties
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.

        Returns
        -------
        Baxis
        """
        super().__init__("baxis")
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
The first argument to the plotly.graph_objs.carpet.Baxis
constructor must be a dict or
an instance of :class:`plotly.graph_objs.carpet.Baxis`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("arraydtick", arg, arraydtick)
        self._set_property("arraytick0", arg, arraytick0)
        self._set_property("autorange", arg, autorange)
        self._set_property("autotypenumbers", arg, autotypenumbers)
        self._set_property("categoryarray", arg, categoryarray)
        self._set_property("categoryarraysrc", arg, categoryarraysrc)
        self._set_property("categoryorder", arg, categoryorder)
        self._set_property("cheatertype", arg, cheatertype)
        self._set_property("color", arg, color)
        self._set_property("dtick", arg, dtick)
        self._set_property("endline", arg, endline)
        self._set_property("endlinecolor", arg, endlinecolor)
        self._set_property("endlinewidth", arg, endlinewidth)
        self._set_property("exponentformat", arg, exponentformat)
        self._set_property("fixedrange", arg, fixedrange)
        self._set_property("gridcolor", arg, gridcolor)
        self._set_property("griddash", arg, griddash)
        self._set_property("gridwidth", arg, gridwidth)
        self._set_property("labelalias", arg, labelalias)
        self._set_property("labelpadding", arg, labelpadding)
        self._set_property("labelprefix", arg, labelprefix)
        self._set_property("labelsuffix", arg, labelsuffix)
        self._set_property("linecolor", arg, linecolor)
        self._set_property("linewidth", arg, linewidth)
        self._set_property("minexponent", arg, minexponent)
        self._set_property("minorgridcolor", arg, minorgridcolor)
        self._set_property("minorgridcount", arg, minorgridcount)
        self._set_property("minorgriddash", arg, minorgriddash)
        self._set_property("minorgridwidth", arg, minorgridwidth)
        self._set_property("nticks", arg, nticks)
        self._set_property("range", arg, range)
        self._set_property("rangemode", arg, rangemode)
        self._set_property("separatethousands", arg, separatethousands)
        self._set_property("showexponent", arg, showexponent)
        self._set_property("showgrid", arg, showgrid)
        self._set_property("showline", arg, showline)
        self._set_property("showticklabels", arg, showticklabels)
        self._set_property("showtickprefix", arg, showtickprefix)
        self._set_property("showticksuffix", arg, showticksuffix)
        self._set_property("smoothing", arg, smoothing)
        self._set_property("startline", arg, startline)
        self._set_property("startlinecolor", arg, startlinecolor)
        self._set_property("startlinewidth", arg, startlinewidth)
        self._set_property("tick0", arg, tick0)
        self._set_property("tickangle", arg, tickangle)
        self._set_property("tickfont", arg, tickfont)
        self._set_property("tickformat", arg, tickformat)
        self._set_property("tickformatstops", arg, tickformatstops)
        self._set_property("tickformatstopdefaults", arg, tickformatstopdefaults)
        self._set_property("tickmode", arg, tickmode)
        self._set_property("tickprefix", arg, tickprefix)
        self._set_property("ticksuffix", arg, ticksuffix)
        self._set_property("ticktext", arg, ticktext)
        self._set_property("ticktextsrc", arg, ticktextsrc)
        self._set_property("tickvals", arg, tickvals)
        self._set_property("tickvalssrc", arg, tickvalssrc)
        self._set_property("title", arg, title)
        self._set_property("type", arg, type)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
