#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):
    _parent_path_str = "box"
    _path_str = "box.marker"
    _valid_props = {
        "angle",
        "color",
        "line",
        "opacity",
        "outliercolor",
        "size",
        "symbol",
    }

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
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.box.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.box.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

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

    @property
    def outliercolor(self):
        """
        Sets the color of the outlier sample points.

        The 'outliercolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["outliercolor"]

    @outliercolor.setter
    def outliercolor(self, val):
        self["outliercolor"] = val

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
            :class:`plotly.graph_objects.box.marker.Line` instance
            or dict with compatible properties
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
            an instance of :class:`plotly.graph_objs.box.Marker`
        angle
            Sets the marker angle in respect to `angleref`.
        color
            Sets the marker color. It accepts either a specific
            color or an array of numbers that are mapped to the
            colorscale relative to the max and min values of the
            array or relative to `marker.cmin` and `marker.cmax` if
            set.
        line
            :class:`plotly.graph_objects.box.marker.Line` instance
            or dict with compatible properties
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
        super().__init__("marker")
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
The first argument to the plotly.graph_objs.box.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.box.Marker`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("angle", arg, angle)
        self._set_property("color", arg, color)
        self._set_property("line", arg, line)
        self._set_property("opacity", arg, opacity)
        self._set_property("outliercolor", arg, outliercolor)
        self._set_property("size", arg, size)
        self._set_property("symbol", arg, symbol)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
