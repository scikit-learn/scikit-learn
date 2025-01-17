from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Aaxis(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "carpet"
    _path_str = "carpet.aaxis"
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
        "titlefont",
        "titleoffset",
        "type",
    }

    # arraydtick
    # ----------
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

    # arraytick0
    # ----------
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

    # autorange
    # ---------
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

    # cheatertype
    # -----------
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

    # dtick
    # -----
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

    # endline
    # -------
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

    # endlinecolor
    # ------------
    @property
    def endlinecolor(self):
        """
        Sets the line color of the end line.

        The 'endlinecolor' property is a color and may be specified as:
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
        return self["endlinecolor"]

    @endlinecolor.setter
    def endlinecolor(self, val):
        self["endlinecolor"] = val

    # endlinewidth
    # ------------
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
        Sets the axis line color.

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

    # labelpadding
    # ------------
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

    # labelprefix
    # -----------
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

    # labelsuffix
    # -----------
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

    # minexponent
    # -----------
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

    # minorgridcolor
    # --------------
    @property
    def minorgridcolor(self):
        """
        Sets the color of the grid lines.

        The 'minorgridcolor' property is a color and may be specified as:
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
        return self["minorgridcolor"]

    @minorgridcolor.setter
    def minorgridcolor(self, val):
        self["minorgridcolor"] = val

    # minorgridcount
    # --------------
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

    # minorgriddash
    # -------------
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

    # minorgridwidth
    # --------------
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

    # rangemode
    # ---------
    @property
    def rangemode(self):
        """
        If "normal", the range is computed in relation to the extrema
        of the input data. If *tozero*`, the range extends to 0,
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

    # showticklabels
    # --------------
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

    # smoothing
    # ---------
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

    # startline
    # ---------
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

    # startlinecolor
    # --------------
    @property
    def startlinecolor(self):
        """
        Sets the line color of the start line.

        The 'startlinecolor' property is a color and may be specified as:
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
        return self["startlinecolor"]

    @startlinecolor.setter
    def startlinecolor(self, val):
        self["startlinecolor"] = val

    # startlinewidth
    # --------------
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

    # tick0
    # -----
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

    # tickfont
    # --------
    @property
    def tickfont(self):
        """
        Sets the tick font.

        The 'tickfont' property is an instance of Tickfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.aaxis.Tickfont`
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
        plotly.graph_objs.carpet.aaxis.Tickfont
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
          - A list or tuple of instances of plotly.graph_objs.carpet.aaxis.Tickformatstop
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
        tuple[plotly.graph_objs.carpet.aaxis.Tickformatstop]
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
        layout.template.data.carpet.aaxis.tickformatstopdefaults), sets
        the default property values to use for elements of
        carpet.aaxis.tickformatstops

        The 'tickformatstopdefaults' property is an instance of Tickformatstop
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.aaxis.Tickformatstop`
          - A dict of string/value properties that will be passed
            to the Tickformatstop constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.carpet.aaxis.Tickformatstop
        """
        return self["tickformatstopdefaults"]

    @tickformatstopdefaults.setter
    def tickformatstopdefaults(self, val):
        self["tickformatstopdefaults"] = val

    # tickmode
    # --------
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

    # title
    # -----
    @property
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.aaxis.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

            Supported dict properties:

                font
                    Sets this axis' title font. Note that the
                    title's font used to be set by the now
                    deprecated `titlefont` attribute.
                offset
                    An additional amount by which to offset the
                    title from the tick labels, given in pixels.
                    Note that this used to be set by the now
                    deprecated `titleoffset` attribute.
                text
                    Sets the title of this axis. Note that before
                    the existence of `title.text`, the title's
                    contents used to be defined as the `title`
                    attribute itself. This behavior has been
                    deprecated.

        Returns
        -------
        plotly.graph_objs.carpet.aaxis.Title
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
        Deprecated: Please use carpet.aaxis.title.font instead. Sets
        this axis' title font. Note that the title's font used to be
        set by the now deprecated `titlefont` attribute.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.aaxis.title.Font`
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

    # titleoffset
    # -----------
    @property
    def titleoffset(self):
        """
        Deprecated: Please use carpet.aaxis.title.offset instead. An
        additional amount by which to offset the title from the tick
        labels, given in pixels. Note that this used to be set by the
        now deprecated `titleoffset` attribute.

        The 'offset' property is a number and may be specified as:
          - An int or float

        Returns
        -------

        """
        return self["titleoffset"]

    @titleoffset.setter
    def titleoffset(self, val):
        self["titleoffset"] = val

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
                ['-', 'linear', 'date', 'category']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    # Self properties description
    # ---------------------------
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
            extrema of the input data. If *tozero*`, the range
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
            A tuple of :class:`plotly.graph_objects.carpet.aaxis.Ti
            ckformatstop` instances or dicts with compatible
            properties
        tickformatstopdefaults
            When used in a template (as layout.template.data.carpet
            .aaxis.tickformatstopdefaults), sets the default
            property values to use for elements of
            carpet.aaxis.tickformatstops
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
            :class:`plotly.graph_objects.carpet.aaxis.Title`
            instance or dict with compatible properties
        titlefont
            Deprecated: Please use carpet.aaxis.title.font instead.
            Sets this axis' title font. Note that the title's font
            used to be set by the now deprecated `titlefont`
            attribute.
        titleoffset
            Deprecated: Please use carpet.aaxis.title.offset
            instead. An additional amount by which to offset the
            title from the tick labels, given in pixels. Note that
            this used to be set by the now deprecated `titleoffset`
            attribute.
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.
        """

    _mapped_properties = {
        "titlefont": ("title", "font"),
        "titleoffset": ("title", "offset"),
    }

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
        titlefont=None,
        titleoffset=None,
        type=None,
        **kwargs,
    ):
        """
        Construct a new Aaxis object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.carpet.Aaxis`
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
            extrema of the input data. If *tozero*`, the range
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
            A tuple of :class:`plotly.graph_objects.carpet.aaxis.Ti
            ckformatstop` instances or dicts with compatible
            properties
        tickformatstopdefaults
            When used in a template (as layout.template.data.carpet
            .aaxis.tickformatstopdefaults), sets the default
            property values to use for elements of
            carpet.aaxis.tickformatstops
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
            :class:`plotly.graph_objects.carpet.aaxis.Title`
            instance or dict with compatible properties
        titlefont
            Deprecated: Please use carpet.aaxis.title.font instead.
            Sets this axis' title font. Note that the title's font
            used to be set by the now deprecated `titlefont`
            attribute.
        titleoffset
            Deprecated: Please use carpet.aaxis.title.offset
            instead. An additional amount by which to offset the
            title from the tick labels, given in pixels. Note that
            this used to be set by the now deprecated `titleoffset`
            attribute.
        type
            Sets the axis type. By default, plotly attempts to
            determined the axis type by looking into the data of
            the traces that referenced the axis in question.

        Returns
        -------
        Aaxis
        """
        super(Aaxis, self).__init__("aaxis")

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
The first argument to the plotly.graph_objs.carpet.Aaxis
constructor must be a dict or
an instance of :class:`plotly.graph_objs.carpet.Aaxis`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("arraydtick", None)
        _v = arraydtick if arraydtick is not None else _v
        if _v is not None:
            self["arraydtick"] = _v
        _v = arg.pop("arraytick0", None)
        _v = arraytick0 if arraytick0 is not None else _v
        if _v is not None:
            self["arraytick0"] = _v
        _v = arg.pop("autorange", None)
        _v = autorange if autorange is not None else _v
        if _v is not None:
            self["autorange"] = _v
        _v = arg.pop("autotypenumbers", None)
        _v = autotypenumbers if autotypenumbers is not None else _v
        if _v is not None:
            self["autotypenumbers"] = _v
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
        _v = arg.pop("cheatertype", None)
        _v = cheatertype if cheatertype is not None else _v
        if _v is not None:
            self["cheatertype"] = _v
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("dtick", None)
        _v = dtick if dtick is not None else _v
        if _v is not None:
            self["dtick"] = _v
        _v = arg.pop("endline", None)
        _v = endline if endline is not None else _v
        if _v is not None:
            self["endline"] = _v
        _v = arg.pop("endlinecolor", None)
        _v = endlinecolor if endlinecolor is not None else _v
        if _v is not None:
            self["endlinecolor"] = _v
        _v = arg.pop("endlinewidth", None)
        _v = endlinewidth if endlinewidth is not None else _v
        if _v is not None:
            self["endlinewidth"] = _v
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
        _v = arg.pop("labelalias", None)
        _v = labelalias if labelalias is not None else _v
        if _v is not None:
            self["labelalias"] = _v
        _v = arg.pop("labelpadding", None)
        _v = labelpadding if labelpadding is not None else _v
        if _v is not None:
            self["labelpadding"] = _v
        _v = arg.pop("labelprefix", None)
        _v = labelprefix if labelprefix is not None else _v
        if _v is not None:
            self["labelprefix"] = _v
        _v = arg.pop("labelsuffix", None)
        _v = labelsuffix if labelsuffix is not None else _v
        if _v is not None:
            self["labelsuffix"] = _v
        _v = arg.pop("linecolor", None)
        _v = linecolor if linecolor is not None else _v
        if _v is not None:
            self["linecolor"] = _v
        _v = arg.pop("linewidth", None)
        _v = linewidth if linewidth is not None else _v
        if _v is not None:
            self["linewidth"] = _v
        _v = arg.pop("minexponent", None)
        _v = minexponent if minexponent is not None else _v
        if _v is not None:
            self["minexponent"] = _v
        _v = arg.pop("minorgridcolor", None)
        _v = minorgridcolor if minorgridcolor is not None else _v
        if _v is not None:
            self["minorgridcolor"] = _v
        _v = arg.pop("minorgridcount", None)
        _v = minorgridcount if minorgridcount is not None else _v
        if _v is not None:
            self["minorgridcount"] = _v
        _v = arg.pop("minorgriddash", None)
        _v = minorgriddash if minorgriddash is not None else _v
        if _v is not None:
            self["minorgriddash"] = _v
        _v = arg.pop("minorgridwidth", None)
        _v = minorgridwidth if minorgridwidth is not None else _v
        if _v is not None:
            self["minorgridwidth"] = _v
        _v = arg.pop("nticks", None)
        _v = nticks if nticks is not None else _v
        if _v is not None:
            self["nticks"] = _v
        _v = arg.pop("range", None)
        _v = range if range is not None else _v
        if _v is not None:
            self["range"] = _v
        _v = arg.pop("rangemode", None)
        _v = rangemode if rangemode is not None else _v
        if _v is not None:
            self["rangemode"] = _v
        _v = arg.pop("separatethousands", None)
        _v = separatethousands if separatethousands is not None else _v
        if _v is not None:
            self["separatethousands"] = _v
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
        _v = arg.pop("smoothing", None)
        _v = smoothing if smoothing is not None else _v
        if _v is not None:
            self["smoothing"] = _v
        _v = arg.pop("startline", None)
        _v = startline if startline is not None else _v
        if _v is not None:
            self["startline"] = _v
        _v = arg.pop("startlinecolor", None)
        _v = startlinecolor if startlinecolor is not None else _v
        if _v is not None:
            self["startlinecolor"] = _v
        _v = arg.pop("startlinewidth", None)
        _v = startlinewidth if startlinewidth is not None else _v
        if _v is not None:
            self["startlinewidth"] = _v
        _v = arg.pop("tick0", None)
        _v = tick0 if tick0 is not None else _v
        if _v is not None:
            self["tick0"] = _v
        _v = arg.pop("tickangle", None)
        _v = tickangle if tickangle is not None else _v
        if _v is not None:
            self["tickangle"] = _v
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
        _v = arg.pop("tickmode", None)
        _v = tickmode if tickmode is not None else _v
        if _v is not None:
            self["tickmode"] = _v
        _v = arg.pop("tickprefix", None)
        _v = tickprefix if tickprefix is not None else _v
        if _v is not None:
            self["tickprefix"] = _v
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
        _v = arg.pop("title", None)
        _v = title if title is not None else _v
        if _v is not None:
            self["title"] = _v
        _v = arg.pop("titlefont", None)
        _v = titlefont if titlefont is not None else _v
        if _v is not None:
            self["titlefont"] = _v
        _v = arg.pop("titleoffset", None)
        _v = titleoffset if titleoffset is not None else _v
        if _v is not None:
            self["titleoffset"] = _v
        _v = arg.pop("type", None)
        _v = type if type is not None else _v
        if _v is not None:
            self["type"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
