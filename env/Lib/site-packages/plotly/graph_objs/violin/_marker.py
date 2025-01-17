from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "violin"
    _path_str = "violin.marker"
    _valid_props = {
        "angle",
        "color",
        "line",
        "opacity",
        "outliercolor",
        "size",
        "symbol",
    }

    # angle
    # -----
    @property
    def angle(self):
        """
        Sets the marker angle in respect to `angleref`.

        The 'angle' property is a angle (in degrees) that may be
        specified as a number between -180 and 180.
        Numeric values outside this range are converted to the equivalent value
        (e.g. 270 is converted to -90).

        Returns
        -------
        int|float
        """
        return self["angle"]

    @angle.setter
    def angle(self, val):
        self["angle"] = val

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

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.violin.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the marker.line color. It accepts either a
                    specific color or an array of numbers that are
                    mapped to the colorscale relative to the max
                    and min values of the array or relative to
                    `marker.line.cmin` and `marker.line.cmax` if
                    set.
                outliercolor
                    Sets the border line color of the outlier
                    sample points. Defaults to marker.color
                outlierwidth
                    Sets the border line width (in px) of the
                    outlier sample points.
                width
                    Sets the width (in px) of the lines bounding
                    the marker points.

        Returns
        -------
        plotly.graph_objs.violin.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # opacity
    # -------
    @property
    def opacity(self):
        """
        Sets the marker opacity.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    # outliercolor
    # ------------
    @property
    def outliercolor(self):
        """
        Sets the color of the outlier sample points.

        The 'outliercolor' property is a color and may be specified as:
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
        return self["outliercolor"]

    @outliercolor.setter
    def outliercolor(self, val):
        self["outliercolor"] = val

    # size
    # ----
    @property
    def size(self):
        """
        Sets the marker size (in px).

        The 'size' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

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

        Returns
        -------
        Any
        """
        return self["symbol"]

    @symbol.setter
    def symbol(self, val):
        self["symbol"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        angle
            Sets the marker angle in respect to `angleref`.
        color
            Sets the marker color. It accepts either a specific
            color or an array of numbers that are mapped to the
            colorscale relative to the max and min values of the
            array or relative to `marker.cmin` and `marker.cmax` if
            set.
        line
            :class:`plotly.graph_objects.violin.marker.Line`
            instance or dict with compatible properties
        opacity
            Sets the marker opacity.
        outliercolor
            Sets the color of the outlier sample points.
        size
            Sets the marker size (in px).
        symbol
            Sets the marker symbol type. Adding 100 is equivalent
            to appending "-open" to a symbol name. Adding 200 is
            equivalent to appending "-dot" to a symbol name. Adding
            300 is equivalent to appending "-open-dot" or "dot-
            open" to a symbol name.
        """

    def __init__(
        self,
        arg=None,
        angle=None,
        color=None,
        line=None,
        opacity=None,
        outliercolor=None,
        size=None,
        symbol=None,
        **kwargs,
    ):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.violin.Marker`
        angle
            Sets the marker angle in respect to `angleref`.
        color
            Sets the marker color. It accepts either a specific
            color or an array of numbers that are mapped to the
            colorscale relative to the max and min values of the
            array or relative to `marker.cmin` and `marker.cmax` if
            set.
        line
            :class:`plotly.graph_objects.violin.marker.Line`
            instance or dict with compatible properties
        opacity
            Sets the marker opacity.
        outliercolor
            Sets the color of the outlier sample points.
        size
            Sets the marker size (in px).
        symbol
            Sets the marker symbol type. Adding 100 is equivalent
            to appending "-open" to a symbol name. Adding 200 is
            equivalent to appending "-dot" to a symbol name. Adding
            300 is equivalent to appending "-open-dot" or "dot-
            open" to a symbol name.

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
The first argument to the plotly.graph_objs.violin.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.violin.Marker`"""
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
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("opacity", None)
        _v = opacity if opacity is not None else _v
        if _v is not None:
            self["opacity"] = _v
        _v = arg.pop("outliercolor", None)
        _v = outliercolor if outliercolor is not None else _v
        if _v is not None:
            self["outliercolor"] = _v
        _v = arg.pop("size", None)
        _v = size if size is not None else _v
        if _v is not None:
            self["size"] = _v
        _v = arg.pop("symbol", None)
        _v = symbol if symbol is not None else _v
        if _v is not None:
            self["symbol"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
