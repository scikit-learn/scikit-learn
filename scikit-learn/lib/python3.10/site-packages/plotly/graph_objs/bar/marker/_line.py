#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Line(_BaseTraceHierarchyType):
    _parent_path_str = "bar.marker"
    _path_str = "bar.marker.line"
    _valid_props = {
        "autocolorscale",
        "cauto",
        "cmax",
        "cmid",
        "cmin",
        "color",
        "coloraxis",
        "colorscale",
        "colorsrc",
        "reversescale",
        "width",
        "widthsrc",
    }

    @property
    def autocolorscale(self):
        """
        Determines whether the colorscale is a default palette
        (`autocolorscale: true`) or the palette determined by
        `marker.line.colorscale`. Has an effect only if in
        `marker.line.color` is set to a numerical array. In case
        `colorscale` is unspecified or `autocolorscale` is true, the
        default palette will be chosen according to whether numbers in
        the `color` array are all positive, all negative or mixed.

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
        respect to the input data (here in `marker.line.color`) or the
        bounds set in `marker.line.cmin` and `marker.line.cmax` Has an
        effect only if in `marker.line.color` is set to a numerical
        array. Defaults to `false` when `marker.line.cmin` and
        `marker.line.cmax` are set by the user.

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
        in `marker.line.color` is set to a numerical array. Value
        should have the same units as in `marker.line.color` and if
        set, `marker.line.cmin` must be set as well.

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
        Sets the mid-point of the color domain by scaling
        `marker.line.cmin` and/or `marker.line.cmax` to be equidistant
        to this point. Has an effect only if in `marker.line.color` is
        set to a numerical array. Value should have the same units as
        in `marker.line.color`. Has no effect when `marker.line.cauto`
        is `false`.

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
        in `marker.line.color` is set to a numerical array. Value
        should have the same units as in `marker.line.color` and if
        set, `marker.line.cmax` must be set as well.

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
        Sets the marker.line color. It accepts either a specific color
        or an array of numbers that are mapped to the colorscale
        relative to the max and min values of the array or relative to
        `marker.line.cmin` and `marker.line.cmax` if set.

        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list
          - A number that will be interpreted as a color
            according to bar.marker.line.colorscale
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
    def colorscale(self):
        """
        Sets the colorscale. Has an effect only if in
        `marker.line.color` is set to a numerical array. The colorscale
        must be an array containing arrays mapping a normalized value
        to an rgb, rgba, hex, hsl, hsv, or named color string. At
        minimum, a mapping for the lowest (0) and highest (1) values
        are required. For example, `[[0, 'rgb(0,0,255)'], [1,
        'rgb(255,0,0)']]`. To control the bounds of the colorscale in
        color space, use `marker.line.cmin` and `marker.line.cmax`.
        Alternatively, `colorscale` may be a palette name string of the
        following list: Blackbody,Bluered,Blues,Cividis,Earth,Electric,
        Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,RdBu,Reds,Viridis,
        YlGnBu,YlOrRd.

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
    def reversescale(self):
        """
        Reverses the color mapping if true. Has an effect only if in
        `marker.line.color` is set to a numerical array. If true,
        `marker.line.cmin` will correspond to the last color in the
        array and `marker.line.cmax` will correspond to the first
        color.

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
    def width(self):
        """
        Sets the width (in px) of the lines bounding the marker points.

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [0, inf]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    @property
    def widthsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `width`.

        The 'widthsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["widthsrc"]

    @widthsrc.setter
    def widthsrc(self, val):
        self["widthsrc"] = val

    @property
    def _prop_descriptions(self):
        return """\
        autocolorscale
            Determines whether the colorscale is a default palette
            (`autocolorscale: true`) or the palette determined by
            `marker.line.colorscale`. Has an effect only if in
            `marker.line.color` is set to a numerical array. In
            case `colorscale` is unspecified or `autocolorscale` is
            true, the default palette will be chosen according to
            whether numbers in the `color` array are all positive,
            all negative or mixed.
        cauto
            Determines whether or not the color domain is computed
            with respect to the input data (here in
            `marker.line.color`) or the bounds set in
            `marker.line.cmin` and `marker.line.cmax` Has an effect
            only if in `marker.line.color` is set to a numerical
            array. Defaults to `false` when `marker.line.cmin` and
            `marker.line.cmax` are set by the user.
        cmax
            Sets the upper bound of the color domain. Has an effect
            only if in `marker.line.color` is set to a numerical
            array. Value should have the same units as in
            `marker.line.color` and if set, `marker.line.cmin` must
            be set as well.
        cmid
            Sets the mid-point of the color domain by scaling
            `marker.line.cmin` and/or `marker.line.cmax` to be
            equidistant to this point. Has an effect only if in
            `marker.line.color` is set to a numerical array. Value
            should have the same units as in `marker.line.color`.
            Has no effect when `marker.line.cauto` is `false`.
        cmin
            Sets the lower bound of the color domain. Has an effect
            only if in `marker.line.color` is set to a numerical
            array. Value should have the same units as in
            `marker.line.color` and if set, `marker.line.cmax` must
            be set as well.
        color
            Sets the marker.line color. It accepts either a
            specific color or an array of numbers that are mapped
            to the colorscale relative to the max and min values of
            the array or relative to `marker.line.cmin` and
            `marker.line.cmax` if set.
        coloraxis
            Sets a reference to a shared color axis. References to
            these shared color axes are "coloraxis", "coloraxis2",
            "coloraxis3", etc. Settings for these shared color axes
            are set in the layout, under `layout.coloraxis`,
            `layout.coloraxis2`, etc. Note that multiple color
            scales can be linked to the same color axis.
        colorscale
            Sets the colorscale. Has an effect only if in
            `marker.line.color` is set to a numerical array. The
            colorscale must be an array containing arrays mapping a
            normalized value to an rgb, rgba, hex, hsl, hsv, or
            named color string. At minimum, a mapping for the
            lowest (0) and highest (1) values are required. For
            example, `[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`.
            To control the bounds of the colorscale in color space,
            use `marker.line.cmin` and `marker.line.cmax`.
            Alternatively, `colorscale` may be a palette name
            string of the following list: Blackbody,Bluered,Blues,C
            ividis,Earth,Electric,Greens,Greys,Hot,Jet,Picnic,Portl
            and,Rainbow,RdBu,Reds,Viridis,YlGnBu,YlOrRd.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        reversescale
            Reverses the color mapping if true. Has an effect only
            if in `marker.line.color` is set to a numerical array.
            If true, `marker.line.cmin` will correspond to the last
            color in the array and `marker.line.cmax` will
            correspond to the first color.
        width
            Sets the width (in px) of the lines bounding the marker
            points.
        widthsrc
            Sets the source reference on Chart Studio Cloud for
            `width`.
        """

    def __init__(
        self,
        arg=None,
        autocolorscale=None,
        cauto=None,
        cmax=None,
        cmid=None,
        cmin=None,
        color=None,
        coloraxis=None,
        colorscale=None,
        colorsrc=None,
        reversescale=None,
        width=None,
        widthsrc=None,
        **kwargs,
    ):
        """
        Construct a new Line object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.bar.marker.Line`
        autocolorscale
            Determines whether the colorscale is a default palette
            (`autocolorscale: true`) or the palette determined by
            `marker.line.colorscale`. Has an effect only if in
            `marker.line.color` is set to a numerical array. In
            case `colorscale` is unspecified or `autocolorscale` is
            true, the default palette will be chosen according to
            whether numbers in the `color` array are all positive,
            all negative or mixed.
        cauto
            Determines whether or not the color domain is computed
            with respect to the input data (here in
            `marker.line.color`) or the bounds set in
            `marker.line.cmin` and `marker.line.cmax` Has an effect
            only if in `marker.line.color` is set to a numerical
            array. Defaults to `false` when `marker.line.cmin` and
            `marker.line.cmax` are set by the user.
        cmax
            Sets the upper bound of the color domain. Has an effect
            only if in `marker.line.color` is set to a numerical
            array. Value should have the same units as in
            `marker.line.color` and if set, `marker.line.cmin` must
            be set as well.
        cmid
            Sets the mid-point of the color domain by scaling
            `marker.line.cmin` and/or `marker.line.cmax` to be
            equidistant to this point. Has an effect only if in
            `marker.line.color` is set to a numerical array. Value
            should have the same units as in `marker.line.color`.
            Has no effect when `marker.line.cauto` is `false`.
        cmin
            Sets the lower bound of the color domain. Has an effect
            only if in `marker.line.color` is set to a numerical
            array. Value should have the same units as in
            `marker.line.color` and if set, `marker.line.cmax` must
            be set as well.
        color
            Sets the marker.line color. It accepts either a
            specific color or an array of numbers that are mapped
            to the colorscale relative to the max and min values of
            the array or relative to `marker.line.cmin` and
            `marker.line.cmax` if set.
        coloraxis
            Sets a reference to a shared color axis. References to
            these shared color axes are "coloraxis", "coloraxis2",
            "coloraxis3", etc. Settings for these shared color axes
            are set in the layout, under `layout.coloraxis`,
            `layout.coloraxis2`, etc. Note that multiple color
            scales can be linked to the same color axis.
        colorscale
            Sets the colorscale. Has an effect only if in
            `marker.line.color` is set to a numerical array. The
            colorscale must be an array containing arrays mapping a
            normalized value to an rgb, rgba, hex, hsl, hsv, or
            named color string. At minimum, a mapping for the
            lowest (0) and highest (1) values are required. For
            example, `[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`.
            To control the bounds of the colorscale in color space,
            use `marker.line.cmin` and `marker.line.cmax`.
            Alternatively, `colorscale` may be a palette name
            string of the following list: Blackbody,Bluered,Blues,C
            ividis,Earth,Electric,Greens,Greys,Hot,Jet,Picnic,Portl
            and,Rainbow,RdBu,Reds,Viridis,YlGnBu,YlOrRd.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        reversescale
            Reverses the color mapping if true. Has an effect only
            if in `marker.line.color` is set to a numerical array.
            If true, `marker.line.cmin` will correspond to the last
            color in the array and `marker.line.cmax` will
            correspond to the first color.
        width
            Sets the width (in px) of the lines bounding the marker
            points.
        widthsrc
            Sets the source reference on Chart Studio Cloud for
            `width`.

        Returns
        -------
        Line
        """
        super().__init__("line")
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
The first argument to the plotly.graph_objs.bar.marker.Line
constructor must be a dict or
an instance of :class:`plotly.graph_objs.bar.marker.Line`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("autocolorscale", arg, autocolorscale)
        self._set_property("cauto", arg, cauto)
        self._set_property("cmax", arg, cmax)
        self._set_property("cmid", arg, cmid)
        self._set_property("cmin", arg, cmin)
        self._set_property("color", arg, color)
        self._set_property("coloraxis", arg, coloraxis)
        self._set_property("colorscale", arg, colorscale)
        self._set_property("colorsrc", arg, colorsrc)
        self._set_property("reversescale", arg, reversescale)
        self._set_property("width", arg, width)
        self._set_property("widthsrc", arg, widthsrc)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
