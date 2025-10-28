#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):
    _parent_path_str = "scatterpolar"
    _path_str = "scatterpolar.marker"
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
          - A number that will be interpreted as a color
            according to scatterpolar.marker.colorscale
          - A list or array of any of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

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

    @property
    def colorbar(self):
        """
        The 'colorbar' property is an instance of ColorBar
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatterpolar.marker.ColorBar`
          - A dict of string/value properties that will be passed
            to the ColorBar constructor

        Returns
        -------
        plotly.graph_objs.scatterpolar.marker.ColorBar
        """
        return self["colorbar"]

    @colorbar.setter
    def colorbar(self, val):
        self["colorbar"] = val

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

    @property
    def gradient(self):
        """
        The 'gradient' property is an instance of Gradient
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatterpolar.marker.Gradient`
          - A dict of string/value properties that will be passed
            to the Gradient constructor

        Returns
        -------
        plotly.graph_objs.scatterpolar.marker.Gradient
        """
        return self["gradient"]

    @gradient.setter
    def gradient(self, val):
        self["gradient"] = val

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatterpolar.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.scatterpolar.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

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
            :class:`plotly.graph_objects.scatterpolar.marker.ColorB
            ar` instance or dict with compatible properties
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
            :class:`plotly.graph_objects.scatterpolar.marker.Gradie
            nt` instance or dict with compatible properties
        line
            :class:`plotly.graph_objects.scatterpolar.marker.Line`
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
            :class:`plotly.graph_objs.scatterpolar.Marker`
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
            :class:`plotly.graph_objects.scatterpolar.marker.ColorB
            ar` instance or dict with compatible properties
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
            :class:`plotly.graph_objects.scatterpolar.marker.Gradie
            nt` instance or dict with compatible properties
        line
            :class:`plotly.graph_objects.scatterpolar.marker.Line`
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
The first argument to the plotly.graph_objs.scatterpolar.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scatterpolar.Marker`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("angle", arg, angle)
        self._set_property("angleref", arg, angleref)
        self._set_property("anglesrc", arg, anglesrc)
        self._set_property("autocolorscale", arg, autocolorscale)
        self._set_property("cauto", arg, cauto)
        self._set_property("cmax", arg, cmax)
        self._set_property("cmid", arg, cmid)
        self._set_property("cmin", arg, cmin)
        self._set_property("color", arg, color)
        self._set_property("coloraxis", arg, coloraxis)
        self._set_property("colorbar", arg, colorbar)
        self._set_property("colorscale", arg, colorscale)
        self._set_property("colorsrc", arg, colorsrc)
        self._set_property("gradient", arg, gradient)
        self._set_property("line", arg, line)
        self._set_property("maxdisplayed", arg, maxdisplayed)
        self._set_property("opacity", arg, opacity)
        self._set_property("opacitysrc", arg, opacitysrc)
        self._set_property("reversescale", arg, reversescale)
        self._set_property("showscale", arg, showscale)
        self._set_property("size", arg, size)
        self._set_property("sizemin", arg, sizemin)
        self._set_property("sizemode", arg, sizemode)
        self._set_property("sizeref", arg, sizeref)
        self._set_property("sizesrc", arg, sizesrc)
        self._set_property("standoff", arg, standoff)
        self._set_property("standoffsrc", arg, standoffsrc)
        self._set_property("symbol", arg, symbol)
        self._set_property("symbolsrc", arg, symbolsrc)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
