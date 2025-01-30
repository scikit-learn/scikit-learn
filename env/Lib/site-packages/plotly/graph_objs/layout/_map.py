from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Map(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.map"
    _valid_props = {
        "bearing",
        "bounds",
        "center",
        "domain",
        "layerdefaults",
        "layers",
        "pitch",
        "style",
        "uirevision",
        "zoom",
    }

    # bearing
    # -------
    @property
    def bearing(self):
        """
        Sets the bearing angle of the map in degrees counter-clockwise
        from North (map.bearing).

        The 'bearing' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["bearing"]

    @bearing.setter
    def bearing(self, val):
        self["bearing"] = val

    # bounds
    # ------
    @property
    def bounds(self):
        """
        The 'bounds' property is an instance of Bounds
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.Bounds`
          - A dict of string/value properties that will be passed
            to the Bounds constructor

            Supported dict properties:

                east
                    Sets the maximum longitude of the map (in
                    degrees East) if `west`, `south` and `north`
                    are declared.
                north
                    Sets the maximum latitude of the map (in
                    degrees North) if `east`, `west` and `south`
                    are declared.
                south
                    Sets the minimum latitude of the map (in
                    degrees North) if `east`, `west` and `north`
                    are declared.
                west
                    Sets the minimum longitude of the map (in
                    degrees East) if `east`, `south` and `north`
                    are declared.

        Returns
        -------
        plotly.graph_objs.layout.map.Bounds
        """
        return self["bounds"]

    @bounds.setter
    def bounds(self, val):
        self["bounds"] = val

    # center
    # ------
    @property
    def center(self):
        """
        The 'center' property is an instance of Center
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.Center`
          - A dict of string/value properties that will be passed
            to the Center constructor

            Supported dict properties:

                lat
                    Sets the latitude of the center of the map (in
                    degrees North).
                lon
                    Sets the longitude of the center of the map (in
                    degrees East).

        Returns
        -------
        plotly.graph_objs.layout.map.Center
        """
        return self["center"]

    @center.setter
    def center(self, val):
        self["center"] = val

    # domain
    # ------
    @property
    def domain(self):
        """
        The 'domain' property is an instance of Domain
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

            Supported dict properties:

                column
                    If there is a layout grid, use the domain for
                    this column in the grid for this map subplot .
                row
                    If there is a layout grid, use the domain for
                    this row in the grid for this map subplot .
                x
                    Sets the horizontal domain of this map subplot
                    (in plot fraction).
                y
                    Sets the vertical domain of this map subplot
                    (in plot fraction).

        Returns
        -------
        plotly.graph_objs.layout.map.Domain
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    # layers
    # ------
    @property
    def layers(self):
        """
        The 'layers' property is a tuple of instances of
        Layer that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.map.Layer
          - A list or tuple of dicts of string/value properties that
            will be passed to the Layer constructor

            Supported dict properties:

                below
                    Determines if the layer will be inserted before
                    the layer with the specified ID. If omitted or
                    set to '', the layer will be inserted above
                    every existing layer.
                circle
                    :class:`plotly.graph_objects.layout.map.layer.C
                    ircle` instance or dict with compatible
                    properties
                color
                    Sets the primary layer color. If `type` is
                    "circle", color corresponds to the circle color
                    (map.layer.paint.circle-color) If `type` is
                    "line", color corresponds to the line color
                    (map.layer.paint.line-color) If `type` is
                    "fill", color corresponds to the fill color
                    (map.layer.paint.fill-color) If `type` is
                    "symbol", color corresponds to the icon color
                    (map.layer.paint.icon-color)
                coordinates
                    Sets the coordinates array contains [longitude,
                    latitude] pairs for the image corners listed in
                    clockwise order: top left, top right, bottom
                    right, bottom left. Only has an effect for
                    "image" `sourcetype`.
                fill
                    :class:`plotly.graph_objects.layout.map.layer.F
                    ill` instance or dict with compatible
                    properties
                line
                    :class:`plotly.graph_objects.layout.map.layer.L
                    ine` instance or dict with compatible
                    properties
                maxzoom
                    Sets the maximum zoom level
                    (map.layer.maxzoom). At zoom levels equal to or
                    greater than the maxzoom, the layer will be
                    hidden.
                minzoom
                    Sets the minimum zoom level
                    (map.layer.minzoom). At zoom levels less than
                    the minzoom, the layer will be hidden.
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
                opacity
                    Sets the opacity of the layer. If `type` is
                    "circle", opacity corresponds to the circle
                    opacity (map.layer.paint.circle-opacity) If
                    `type` is "line", opacity corresponds to the
                    line opacity (map.layer.paint.line-opacity) If
                    `type` is "fill", opacity corresponds to the
                    fill opacity (map.layer.paint.fill-opacity) If
                    `type` is "symbol", opacity corresponds to the
                    icon/text opacity (map.layer.paint.text-
                    opacity)
                source
                    Sets the source data for this layer
                    (map.layer.source). When `sourcetype` is set to
                    "geojson", `source` can be a URL to a GeoJSON
                    or a GeoJSON object. When `sourcetype` is set
                    to "vector" or "raster", `source` can be a URL
                    or an array of tile URLs. When `sourcetype` is
                    set to "image", `source` can be a URL to an
                    image.
                sourceattribution
                    Sets the attribution for this source.
                sourcelayer
                    Specifies the layer to use from a vector tile
                    source (map.layer.source-layer). Required for
                    "vector" source type that supports multiple
                    layers.
                sourcetype
                    Sets the source type for this layer, that is
                    the type of the layer data.
                symbol
                    :class:`plotly.graph_objects.layout.map.layer.S
                    ymbol` instance or dict with compatible
                    properties
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
                type
                    Sets the layer type, that is the how the layer
                    data set in `source` will be rendered With
                    `sourcetype` set to "geojson", the following
                    values are allowed: "circle", "line", "fill"
                    and "symbol". but note that "line" and "fill"
                    are not compatible with Point GeoJSON
                    geometries. With `sourcetype` set to "vector",
                    the following values are allowed:  "circle",
                    "line", "fill" and "symbol". With `sourcetype`
                    set to "raster" or `*image*`, only the "raster"
                    value is allowed.
                visible
                    Determines whether this layer is displayed

        Returns
        -------
        tuple[plotly.graph_objs.layout.map.Layer]
        """
        return self["layers"]

    @layers.setter
    def layers(self, val):
        self["layers"] = val

    # layerdefaults
    # -------------
    @property
    def layerdefaults(self):
        """
        When used in a template (as
        layout.template.layout.map.layerdefaults), sets the default
        property values to use for elements of layout.map.layers

        The 'layerdefaults' property is an instance of Layer
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.Layer`
          - A dict of string/value properties that will be passed
            to the Layer constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.map.Layer
        """
        return self["layerdefaults"]

    @layerdefaults.setter
    def layerdefaults(self, val):
        self["layerdefaults"] = val

    # pitch
    # -----
    @property
    def pitch(self):
        """
        Sets the pitch angle of the map (in degrees, where 0 means
        perpendicular to the surface of the map) (map.pitch).

        The 'pitch' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["pitch"]

    @pitch.setter
    def pitch(self, val):
        self["pitch"] = val

    # style
    # -----
    @property
    def style(self):
        """
        Defines the map layers that are rendered by default below the
        trace layers defined in `data`, which are themselves by default
        rendered below the layers defined in `layout.map.layers`.
        These layers can be defined either explicitly as a Map Style
        object which can contain multiple layer definitions that load
        data from any public or private Tile Map Service (TMS or XYZ)
        or Web Map Service (WMS) or implicitly by using one of the
        built-in style objects which use WMSes or by using a custom
        style URL  Map Style objects are of the form described in the
        MapLibre GL JS documentation available at
        https://maplibre.org/maplibre-style-spec/  The built-in
        plotly.js styles objects are: basic, carto-darkmatter, carto-
        darkmatter-nolabels, carto-positron, carto-positron-nolabels,
        carto-voyager, carto-voyager-nolabels, dark, light, open-
        street-map, outdoors, satellite, satellite-streets, streets,
        white-bg.

        The 'style' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["style"]

    @style.setter
    def style(self, val):
        self["style"] = val

    # uirevision
    # ----------
    @property
    def uirevision(self):
        """
        Controls persistence of user-driven changes in the view:
        `center`, `zoom`, `bearing`, `pitch`. Defaults to
        `layout.uirevision`.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    # zoom
    # ----
    @property
    def zoom(self):
        """
        Sets the zoom level of the map (map.zoom).

        The 'zoom' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["zoom"]

    @zoom.setter
    def zoom(self, val):
        self["zoom"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        bearing
            Sets the bearing angle of the map in degrees counter-
            clockwise from North (map.bearing).
        bounds
            :class:`plotly.graph_objects.layout.map.Bounds`
            instance or dict with compatible properties
        center
            :class:`plotly.graph_objects.layout.map.Center`
            instance or dict with compatible properties
        domain
            :class:`plotly.graph_objects.layout.map.Domain`
            instance or dict with compatible properties
        layers
            A tuple of
            :class:`plotly.graph_objects.layout.map.Layer`
            instances or dicts with compatible properties
        layerdefaults
            When used in a template (as
            layout.template.layout.map.layerdefaults), sets the
            default property values to use for elements of
            layout.map.layers
        pitch
            Sets the pitch angle of the map (in degrees, where 0
            means perpendicular to the surface of the map)
            (map.pitch).
        style
            Defines the map layers that are rendered by default
            below the trace layers defined in `data`, which are
            themselves by default rendered below the layers defined
            in `layout.map.layers`.  These layers can be defined
            either explicitly as a Map Style object which can
            contain multiple layer definitions that load data from
            any public or private Tile Map Service (TMS or XYZ) or
            Web Map Service (WMS) or implicitly by using one of the
            built-in style objects which use WMSes or by using a
            custom style URL  Map Style objects are of the form
            described in the MapLibre GL JS documentation available
            at https://maplibre.org/maplibre-style-spec/  The
            built-in plotly.js styles objects are: basic, carto-
            darkmatter, carto-darkmatter-nolabels, carto-positron,
            carto-positron-nolabels, carto-voyager, carto-voyager-
            nolabels, dark, light, open-street-map, outdoors,
            satellite, satellite-streets, streets, white-bg.
        uirevision
            Controls persistence of user-driven changes in the
            view: `center`, `zoom`, `bearing`, `pitch`. Defaults to
            `layout.uirevision`.
        zoom
            Sets the zoom level of the map (map.zoom).
        """

    def __init__(
        self,
        arg=None,
        bearing=None,
        bounds=None,
        center=None,
        domain=None,
        layers=None,
        layerdefaults=None,
        pitch=None,
        style=None,
        uirevision=None,
        zoom=None,
        **kwargs,
    ):
        """
        Construct a new Map object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Map`
        bearing
            Sets the bearing angle of the map in degrees counter-
            clockwise from North (map.bearing).
        bounds
            :class:`plotly.graph_objects.layout.map.Bounds`
            instance or dict with compatible properties
        center
            :class:`plotly.graph_objects.layout.map.Center`
            instance or dict with compatible properties
        domain
            :class:`plotly.graph_objects.layout.map.Domain`
            instance or dict with compatible properties
        layers
            A tuple of
            :class:`plotly.graph_objects.layout.map.Layer`
            instances or dicts with compatible properties
        layerdefaults
            When used in a template (as
            layout.template.layout.map.layerdefaults), sets the
            default property values to use for elements of
            layout.map.layers
        pitch
            Sets the pitch angle of the map (in degrees, where 0
            means perpendicular to the surface of the map)
            (map.pitch).
        style
            Defines the map layers that are rendered by default
            below the trace layers defined in `data`, which are
            themselves by default rendered below the layers defined
            in `layout.map.layers`.  These layers can be defined
            either explicitly as a Map Style object which can
            contain multiple layer definitions that load data from
            any public or private Tile Map Service (TMS or XYZ) or
            Web Map Service (WMS) or implicitly by using one of the
            built-in style objects which use WMSes or by using a
            custom style URL  Map Style objects are of the form
            described in the MapLibre GL JS documentation available
            at https://maplibre.org/maplibre-style-spec/  The
            built-in plotly.js styles objects are: basic, carto-
            darkmatter, carto-darkmatter-nolabels, carto-positron,
            carto-positron-nolabels, carto-voyager, carto-voyager-
            nolabels, dark, light, open-street-map, outdoors,
            satellite, satellite-streets, streets, white-bg.
        uirevision
            Controls persistence of user-driven changes in the
            view: `center`, `zoom`, `bearing`, `pitch`. Defaults to
            `layout.uirevision`.
        zoom
            Sets the zoom level of the map (map.zoom).

        Returns
        -------
        Map
        """
        super(Map, self).__init__("map")

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
The first argument to the plotly.graph_objs.layout.Map
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Map`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("bearing", None)
        _v = bearing if bearing is not None else _v
        if _v is not None:
            self["bearing"] = _v
        _v = arg.pop("bounds", None)
        _v = bounds if bounds is not None else _v
        if _v is not None:
            self["bounds"] = _v
        _v = arg.pop("center", None)
        _v = center if center is not None else _v
        if _v is not None:
            self["center"] = _v
        _v = arg.pop("domain", None)
        _v = domain if domain is not None else _v
        if _v is not None:
            self["domain"] = _v
        _v = arg.pop("layers", None)
        _v = layers if layers is not None else _v
        if _v is not None:
            self["layers"] = _v
        _v = arg.pop("layerdefaults", None)
        _v = layerdefaults if layerdefaults is not None else _v
        if _v is not None:
            self["layerdefaults"] = _v
        _v = arg.pop("pitch", None)
        _v = pitch if pitch is not None else _v
        if _v is not None:
            self["pitch"] = _v
        _v = arg.pop("style", None)
        _v = style if style is not None else _v
        if _v is not None:
            self["style"] = _v
        _v = arg.pop("uirevision", None)
        _v = uirevision if uirevision is not None else _v
        if _v is not None:
            self["uirevision"] = _v
        _v = arg.pop("zoom", None)
        _v = zoom if zoom is not None else _v
        if _v is not None:
            self["zoom"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
