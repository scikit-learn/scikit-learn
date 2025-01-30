from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Layer(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.map"
    _path_str = "layout.map.layer"
    _valid_props = {
        "below",
        "circle",
        "color",
        "coordinates",
        "fill",
        "line",
        "maxzoom",
        "minzoom",
        "name",
        "opacity",
        "source",
        "sourceattribution",
        "sourcelayer",
        "sourcetype",
        "symbol",
        "templateitemname",
        "type",
        "visible",
    }

    # below
    # -----
    @property
    def below(self):
        """
        Determines if the layer will be inserted before the layer with
        the specified ID. If omitted or set to '', the layer will be
        inserted above every existing layer.

        The 'below' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["below"]

    @below.setter
    def below(self, val):
        self["below"] = val

    # circle
    # ------
    @property
    def circle(self):
        """
        The 'circle' property is an instance of Circle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.layer.Circle`
          - A dict of string/value properties that will be passed
            to the Circle constructor

            Supported dict properties:

                radius
                    Sets the circle radius (map.layer.paint.circle-
                    radius). Has an effect only when `type` is set
                    to "circle".

        Returns
        -------
        plotly.graph_objs.layout.map.layer.Circle
        """
        return self["circle"]

    @circle.setter
    def circle(self, val):
        self["circle"] = val

    # color
    # -----
    @property
    def color(self):
        """
        Sets the primary layer color. If `type` is "circle", color
        corresponds to the circle color (map.layer.paint.circle-color)
        If `type` is "line", color corresponds to the line color
        (map.layer.paint.line-color) If `type` is "fill", color
        corresponds to the fill color (map.layer.paint.fill-color) If
        `type` is "symbol", color corresponds to the icon color
        (map.layer.paint.icon-color)

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

    # coordinates
    # -----------
    @property
    def coordinates(self):
        """
        Sets the coordinates array contains [longitude, latitude] pairs
        for the image corners listed in clockwise order: top left, top
        right, bottom right, bottom left. Only has an effect for
        "image" `sourcetype`.

        The 'coordinates' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["coordinates"]

    @coordinates.setter
    def coordinates(self, val):
        self["coordinates"] = val

    # fill
    # ----
    @property
    def fill(self):
        """
        The 'fill' property is an instance of Fill
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.layer.Fill`
          - A dict of string/value properties that will be passed
            to the Fill constructor

            Supported dict properties:

                outlinecolor
                    Sets the fill outline color
                    (map.layer.paint.fill-outline-color). Has an
                    effect only when `type` is set to "fill".

        Returns
        -------
        plotly.graph_objs.layout.map.layer.Fill
        """
        return self["fill"]

    @fill.setter
    def fill(self, val):
        self["fill"] = val

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.layer.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                dash
                    Sets the length of dashes and gaps
                    (map.layer.paint.line-dasharray). Has an effect
                    only when `type` is set to "line".
                dashsrc
                    Sets the source reference on Chart Studio Cloud
                    for `dash`.
                width
                    Sets the line width (map.layer.paint.line-
                    width). Has an effect only when `type` is set
                    to "line".

        Returns
        -------
        plotly.graph_objs.layout.map.layer.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # maxzoom
    # -------
    @property
    def maxzoom(self):
        """
        Sets the maximum zoom level (map.layer.maxzoom). At zoom levels
        equal to or greater than the maxzoom, the layer will be hidden.

        The 'maxzoom' property is a number and may be specified as:
          - An int or float in the interval [0, 24]

        Returns
        -------
        int|float
        """
        return self["maxzoom"]

    @maxzoom.setter
    def maxzoom(self, val):
        self["maxzoom"] = val

    # minzoom
    # -------
    @property
    def minzoom(self):
        """
        Sets the minimum zoom level (map.layer.minzoom). At zoom levels
        less than the minzoom, the layer will be hidden.

        The 'minzoom' property is a number and may be specified as:
          - An int or float in the interval [0, 24]

        Returns
        -------
        int|float
        """
        return self["minzoom"]

    @minzoom.setter
    def minzoom(self, val):
        self["minzoom"] = val

    # name
    # ----
    @property
    def name(self):
        """
        When used in a template, named items are created in the output
        figure in addition to any items the figure already has in this
        array. You can modify these items in the output figure by
        making your own item with `templateitemname` matching this
        `name` alongside your modifications (including `visible: false`
        or `enabled: false` to hide it). Has no effect outside of a
        template.

        The 'name' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["name"]

    @name.setter
    def name(self, val):
        self["name"] = val

    # opacity
    # -------
    @property
    def opacity(self):
        """
        Sets the opacity of the layer. If `type` is "circle", opacity
        corresponds to the circle opacity (map.layer.paint.circle-
        opacity) If `type` is "line", opacity corresponds to the line
        opacity (map.layer.paint.line-opacity) If `type` is "fill",
        opacity corresponds to the fill opacity (map.layer.paint.fill-
        opacity) If `type` is "symbol", opacity corresponds to the
        icon/text opacity (map.layer.paint.text-opacity)

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

    # source
    # ------
    @property
    def source(self):
        """
        Sets the source data for this layer (map.layer.source). When
        `sourcetype` is set to "geojson", `source` can be a URL to a
        GeoJSON or a GeoJSON object. When `sourcetype` is set to
        "vector" or "raster", `source` can be a URL or an array of tile
        URLs. When `sourcetype` is set to "image", `source` can be a
        URL to an image.

        The 'source' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["source"]

    @source.setter
    def source(self, val):
        self["source"] = val

    # sourceattribution
    # -----------------
    @property
    def sourceattribution(self):
        """
        Sets the attribution for this source.

        The 'sourceattribution' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["sourceattribution"]

    @sourceattribution.setter
    def sourceattribution(self, val):
        self["sourceattribution"] = val

    # sourcelayer
    # -----------
    @property
    def sourcelayer(self):
        """
        Specifies the layer to use from a vector tile source
        (map.layer.source-layer). Required for "vector" source type
        that supports multiple layers.

        The 'sourcelayer' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["sourcelayer"]

    @sourcelayer.setter
    def sourcelayer(self, val):
        self["sourcelayer"] = val

    # sourcetype
    # ----------
    @property
    def sourcetype(self):
        """
        Sets the source type for this layer, that is the type of the
        layer data.

        The 'sourcetype' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['geojson', 'vector', 'raster', 'image']

        Returns
        -------
        Any
        """
        return self["sourcetype"]

    @sourcetype.setter
    def sourcetype(self, val):
        self["sourcetype"] = val

    # symbol
    # ------
    @property
    def symbol(self):
        """
        The 'symbol' property is an instance of Symbol
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.layer.Symbol`
          - A dict of string/value properties that will be passed
            to the Symbol constructor

            Supported dict properties:

                icon
                    Sets the symbol icon image
                    (map.layer.layout.icon-image). Full list:
                    https://www.map.com/maki-icons/
                iconsize
                    Sets the symbol icon size
                    (map.layer.layout.icon-size). Has an effect
                    only when `type` is set to "symbol".
                placement
                    Sets the symbol and/or text placement
                    (map.layer.layout.symbol-placement). If
                    `placement` is "point", the label is placed
                    where the geometry is located If `placement` is
                    "line", the label is placed along the line of
                    the geometry If `placement` is "line-center",
                    the label is placed on the center of the
                    geometry
                text
                    Sets the symbol text (map.layer.layout.text-
                    field).
                textfont
                    Sets the icon text font
                    (color=map.layer.paint.text-color,
                    size=map.layer.layout.text-size). Has an effect
                    only when `type` is set to "symbol".
                textposition
                    Sets the positions of the `text` elements with
                    respects to the (x,y) coordinates.

        Returns
        -------
        plotly.graph_objs.layout.map.layer.Symbol
        """
        return self["symbol"]

    @symbol.setter
    def symbol(self, val):
        self["symbol"] = val

    # templateitemname
    # ----------------
    @property
    def templateitemname(self):
        """
        Used to refer to a named item in this array in the template.
        Named items from the template will be created even without a
        matching item in the input figure, but you can modify one by
        making an item with `templateitemname` matching its `name`,
        alongside your modifications (including `visible: false` or
        `enabled: false` to hide it). If there is no template or no
        matching item, this item will be hidden unless you explicitly
        show it with `visible: true`.

        The 'templateitemname' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["templateitemname"]

    @templateitemname.setter
    def templateitemname(self, val):
        self["templateitemname"] = val

    # type
    # ----
    @property
    def type(self):
        """
        Sets the layer type, that is the how the layer data set in
        `source` will be rendered With `sourcetype` set to "geojson",
        the following values are allowed: "circle", "line", "fill" and
        "symbol". but note that "line" and "fill" are not compatible
        with Point GeoJSON geometries. With `sourcetype` set to
        "vector", the following values are allowed:  "circle", "line",
        "fill" and "symbol". With `sourcetype` set to "raster" or
        `*image*`, only the "raster" value is allowed.

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['circle', 'line', 'fill', 'symbol', 'raster']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines whether this layer is displayed

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

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        below
            Determines if the layer will be inserted before the
            layer with the specified ID. If omitted or set to '',
            the layer will be inserted above every existing layer.
        circle
            :class:`plotly.graph_objects.layout.map.layer.Circle`
            instance or dict with compatible properties
        color
            Sets the primary layer color. If `type` is "circle",
            color corresponds to the circle color
            (map.layer.paint.circle-color) If `type` is "line",
            color corresponds to the line color
            (map.layer.paint.line-color) If `type` is "fill", color
            corresponds to the fill color (map.layer.paint.fill-
            color) If `type` is "symbol", color corresponds to the
            icon color (map.layer.paint.icon-color)
        coordinates
            Sets the coordinates array contains [longitude,
            latitude] pairs for the image corners listed in
            clockwise order: top left, top right, bottom right,
            bottom left. Only has an effect for "image"
            `sourcetype`.
        fill
            :class:`plotly.graph_objects.layout.map.layer.Fill`
            instance or dict with compatible properties
        line
            :class:`plotly.graph_objects.layout.map.layer.Line`
            instance or dict with compatible properties
        maxzoom
            Sets the maximum zoom level (map.layer.maxzoom). At
            zoom levels equal to or greater than the maxzoom, the
            layer will be hidden.
        minzoom
            Sets the minimum zoom level (map.layer.minzoom). At
            zoom levels less than the minzoom, the layer will be
            hidden.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        opacity
            Sets the opacity of the layer. If `type` is "circle",
            opacity corresponds to the circle opacity
            (map.layer.paint.circle-opacity) If `type` is "line",
            opacity corresponds to the line opacity
            (map.layer.paint.line-opacity) If `type` is "fill",
            opacity corresponds to the fill opacity
            (map.layer.paint.fill-opacity) If `type` is "symbol",
            opacity corresponds to the icon/text opacity
            (map.layer.paint.text-opacity)
        source
            Sets the source data for this layer (map.layer.source).
            When `sourcetype` is set to "geojson", `source` can be
            a URL to a GeoJSON or a GeoJSON object. When
            `sourcetype` is set to "vector" or "raster", `source`
            can be a URL or an array of tile URLs. When
            `sourcetype` is set to "image", `source` can be a URL
            to an image.
        sourceattribution
            Sets the attribution for this source.
        sourcelayer
            Specifies the layer to use from a vector tile source
            (map.layer.source-layer). Required for "vector" source
            type that supports multiple layers.
        sourcetype
            Sets the source type for this layer, that is the type
            of the layer data.
        symbol
            :class:`plotly.graph_objects.layout.map.layer.Symbol`
            instance or dict with compatible properties
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        type
            Sets the layer type, that is the how the layer data set
            in `source` will be rendered With `sourcetype` set to
            "geojson", the following values are allowed: "circle",
            "line", "fill" and "symbol". but note that "line" and
            "fill" are not compatible with Point GeoJSON
            geometries. With `sourcetype` set to "vector", the
            following values are allowed:  "circle", "line", "fill"
            and "symbol". With `sourcetype` set to "raster" or
            `*image*`, only the "raster" value is allowed.
        visible
            Determines whether this layer is displayed
        """

    def __init__(
        self,
        arg=None,
        below=None,
        circle=None,
        color=None,
        coordinates=None,
        fill=None,
        line=None,
        maxzoom=None,
        minzoom=None,
        name=None,
        opacity=None,
        source=None,
        sourceattribution=None,
        sourcelayer=None,
        sourcetype=None,
        symbol=None,
        templateitemname=None,
        type=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Layer object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.map.Layer`
        below
            Determines if the layer will be inserted before the
            layer with the specified ID. If omitted or set to '',
            the layer will be inserted above every existing layer.
        circle
            :class:`plotly.graph_objects.layout.map.layer.Circle`
            instance or dict with compatible properties
        color
            Sets the primary layer color. If `type` is "circle",
            color corresponds to the circle color
            (map.layer.paint.circle-color) If `type` is "line",
            color corresponds to the line color
            (map.layer.paint.line-color) If `type` is "fill", color
            corresponds to the fill color (map.layer.paint.fill-
            color) If `type` is "symbol", color corresponds to the
            icon color (map.layer.paint.icon-color)
        coordinates
            Sets the coordinates array contains [longitude,
            latitude] pairs for the image corners listed in
            clockwise order: top left, top right, bottom right,
            bottom left. Only has an effect for "image"
            `sourcetype`.
        fill
            :class:`plotly.graph_objects.layout.map.layer.Fill`
            instance or dict with compatible properties
        line
            :class:`plotly.graph_objects.layout.map.layer.Line`
            instance or dict with compatible properties
        maxzoom
            Sets the maximum zoom level (map.layer.maxzoom). At
            zoom levels equal to or greater than the maxzoom, the
            layer will be hidden.
        minzoom
            Sets the minimum zoom level (map.layer.minzoom). At
            zoom levels less than the minzoom, the layer will be
            hidden.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        opacity
            Sets the opacity of the layer. If `type` is "circle",
            opacity corresponds to the circle opacity
            (map.layer.paint.circle-opacity) If `type` is "line",
            opacity corresponds to the line opacity
            (map.layer.paint.line-opacity) If `type` is "fill",
            opacity corresponds to the fill opacity
            (map.layer.paint.fill-opacity) If `type` is "symbol",
            opacity corresponds to the icon/text opacity
            (map.layer.paint.text-opacity)
        source
            Sets the source data for this layer (map.layer.source).
            When `sourcetype` is set to "geojson", `source` can be
            a URL to a GeoJSON or a GeoJSON object. When
            `sourcetype` is set to "vector" or "raster", `source`
            can be a URL or an array of tile URLs. When
            `sourcetype` is set to "image", `source` can be a URL
            to an image.
        sourceattribution
            Sets the attribution for this source.
        sourcelayer
            Specifies the layer to use from a vector tile source
            (map.layer.source-layer). Required for "vector" source
            type that supports multiple layers.
        sourcetype
            Sets the source type for this layer, that is the type
            of the layer data.
        symbol
            :class:`plotly.graph_objects.layout.map.layer.Symbol`
            instance or dict with compatible properties
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        type
            Sets the layer type, that is the how the layer data set
            in `source` will be rendered With `sourcetype` set to
            "geojson", the following values are allowed: "circle",
            "line", "fill" and "symbol". but note that "line" and
            "fill" are not compatible with Point GeoJSON
            geometries. With `sourcetype` set to "vector", the
            following values are allowed:  "circle", "line", "fill"
            and "symbol". With `sourcetype` set to "raster" or
            `*image*`, only the "raster" value is allowed.
        visible
            Determines whether this layer is displayed

        Returns
        -------
        Layer
        """
        super(Layer, self).__init__("layers")

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
The first argument to the plotly.graph_objs.layout.map.Layer
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.map.Layer`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("below", None)
        _v = below if below is not None else _v
        if _v is not None:
            self["below"] = _v
        _v = arg.pop("circle", None)
        _v = circle if circle is not None else _v
        if _v is not None:
            self["circle"] = _v
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("coordinates", None)
        _v = coordinates if coordinates is not None else _v
        if _v is not None:
            self["coordinates"] = _v
        _v = arg.pop("fill", None)
        _v = fill if fill is not None else _v
        if _v is not None:
            self["fill"] = _v
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("maxzoom", None)
        _v = maxzoom if maxzoom is not None else _v
        if _v is not None:
            self["maxzoom"] = _v
        _v = arg.pop("minzoom", None)
        _v = minzoom if minzoom is not None else _v
        if _v is not None:
            self["minzoom"] = _v
        _v = arg.pop("name", None)
        _v = name if name is not None else _v
        if _v is not None:
            self["name"] = _v
        _v = arg.pop("opacity", None)
        _v = opacity if opacity is not None else _v
        if _v is not None:
            self["opacity"] = _v
        _v = arg.pop("source", None)
        _v = source if source is not None else _v
        if _v is not None:
            self["source"] = _v
        _v = arg.pop("sourceattribution", None)
        _v = sourceattribution if sourceattribution is not None else _v
        if _v is not None:
            self["sourceattribution"] = _v
        _v = arg.pop("sourcelayer", None)
        _v = sourcelayer if sourcelayer is not None else _v
        if _v is not None:
            self["sourcelayer"] = _v
        _v = arg.pop("sourcetype", None)
        _v = sourcetype if sourcetype is not None else _v
        if _v is not None:
            self["sourcetype"] = _v
        _v = arg.pop("symbol", None)
        _v = symbol if symbol is not None else _v
        if _v is not None:
            self["symbol"] = _v
        _v = arg.pop("templateitemname", None)
        _v = templateitemname if templateitemname is not None else _v
        if _v is not None:
            self["templateitemname"] = _v
        _v = arg.pop("type", None)
        _v = type if type is not None else _v
        if _v is not None:
            self["type"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
