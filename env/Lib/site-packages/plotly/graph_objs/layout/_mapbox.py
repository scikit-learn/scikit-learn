from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Mapbox(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.mapbox"
    _valid_props = {
        "accesstoken",
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

    # accesstoken
    # -----------
    @property
    def accesstoken(self):
        """
        Sets the mapbox access token to be used for this mapbox map.
        Alternatively, the mapbox access token can be set in the
        configuration options under `mapboxAccessToken`. Note that
        accessToken are only required when `style` (e.g with values :
        basic, streets, outdoors, light, dark, satellite, satellite-
        streets ) and/or a layout layer references the Mapbox server.

        The 'accesstoken' property is a string and must be specified as:
          - A non-empty string

        Returns
        -------
        str
        """
        return self["accesstoken"]

    @accesstoken.setter
    def accesstoken(self, val):
        self["accesstoken"] = val

    # bearing
    # -------
    @property
    def bearing(self):
        """
        Sets the bearing angle of the map in degrees counter-clockwise
        from North (mapbox.bearing).

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
          - An instance of :class:`plotly.graph_objs.layout.mapbox.Bounds`
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
        plotly.graph_objs.layout.mapbox.Bounds
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
          - An instance of :class:`plotly.graph_objs.layout.mapbox.Center`
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
        plotly.graph_objs.layout.mapbox.Center
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
          - An instance of :class:`plotly.graph_objs.layout.mapbox.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

            Supported dict properties:

                column
                    If there is a layout grid, use the domain for
                    this column in the grid for this mapbox subplot
                    .
                row
                    If there is a layout grid, use the domain for
                    this row in the grid for this mapbox subplot .
                x
                    Sets the horizontal domain of this mapbox
                    subplot (in plot fraction).
                y
                    Sets the vertical domain of this mapbox subplot
                    (in plot fraction).

        Returns
        -------
        plotly.graph_objs.layout.mapbox.Domain
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
          - A list or tuple of instances of plotly.graph_objs.layout.mapbox.Layer
          - A list or tuple of dicts of string/value properties that
            will be passed to the Layer constructor

            Supported dict properties:

                below
                    Determines if the layer will be inserted before
                    the layer with the specified ID. If omitted or
                    set to '', the layer will be inserted above
                    every existing layer.
                circle
                    :class:`plotly.graph_objects.layout.mapbox.laye
                    r.Circle` instance or dict with compatible
                    properties
                color
                    Sets the primary layer color. If `type` is
                    "circle", color corresponds to the circle color
                    (mapbox.layer.paint.circle-color) If `type` is
                    "line", color corresponds to the line color
                    (mapbox.layer.paint.line-color) If `type` is
                    "fill", color corresponds to the fill color
                    (mapbox.layer.paint.fill-color) If `type` is
                    "symbol", color corresponds to the icon color
                    (mapbox.layer.paint.icon-color)
                coordinates
                    Sets the coordinates array contains [longitude,
                    latitude] pairs for the image corners listed in
                    clockwise order: top left, top right, bottom
                    right, bottom left. Only has an effect for
                    "image" `sourcetype`.
                fill
                    :class:`plotly.graph_objects.layout.mapbox.laye
                    r.Fill` instance or dict with compatible
                    properties
                line
                    :class:`plotly.graph_objects.layout.mapbox.laye
                    r.Line` instance or dict with compatible
                    properties
                maxzoom
                    Sets the maximum zoom level
                    (mapbox.layer.maxzoom). At zoom levels equal to
                    or greater than the maxzoom, the layer will be
                    hidden.
                minzoom
                    Sets the minimum zoom level
                    (mapbox.layer.minzoom). At zoom levels less
                    than the minzoom, the layer will be hidden.
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
                    opacity (mapbox.layer.paint.circle-opacity) If
                    `type` is "line", opacity corresponds to the
                    line opacity (mapbox.layer.paint.line-opacity)
                    If `type` is "fill", opacity corresponds to the
                    fill opacity (mapbox.layer.paint.fill-opacity)
                    If `type` is "symbol", opacity corresponds to
                    the icon/text opacity (mapbox.layer.paint.text-
                    opacity)
                source
                    Sets the source data for this layer
                    (mapbox.layer.source). When `sourcetype` is set
                    to "geojson", `source` can be a URL to a
                    GeoJSON or a GeoJSON object. When `sourcetype`
                    is set to "vector" or "raster", `source` can be
                    a URL or an array of tile URLs. When
                    `sourcetype` is set to "image", `source` can be
                    a URL to an image.
                sourceattribution
                    Sets the attribution for this source.
                sourcelayer
                    Specifies the layer to use from a vector tile
                    source (mapbox.layer.source-layer). Required
                    for "vector" source type that supports multiple
                    layers.
                sourcetype
                    Sets the source type for this layer, that is
                    the type of the layer data.
                symbol
                    :class:`plotly.graph_objects.layout.mapbox.laye
                    r.Symbol` instance or dict with compatible
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
        tuple[plotly.graph_objs.layout.mapbox.Layer]
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
        layout.template.layout.mapbox.layerdefaults), sets the default
        property values to use for elements of layout.mapbox.layers

        The 'layerdefaults' property is an instance of Layer
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.mapbox.Layer`
          - A dict of string/value properties that will be passed
            to the Layer constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.mapbox.Layer
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
        perpendicular to the surface of the map) (mapbox.pitch).

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
        rendered below the layers defined in `layout.mapbox.layers`.
        These layers can be defined either explicitly as a Mapbox Style
        object which can contain multiple layer definitions that load
        data from any public or private Tile Map Service (TMS or XYZ)
        or Web Map Service (WMS) or implicitly by using one of the
        built-in style objects which use WMSes which do not require any
        access tokens, or by using a default Mapbox style or custom
        Mapbox style URL, both of which require a Mapbox access token
        Note that Mapbox access token can be set in the `accesstoken`
        attribute or in the `mapboxAccessToken` config option.  Mapbox
        Style objects are of the form described in the Mapbox GL JS
        documentation available at https://docs.mapbox.com/mapbox-gl-
        js/style-spec  The built-in plotly.js styles objects are:
        carto-darkmatter, carto-positron, open-street-map, stamen-
        terrain, stamen-toner, stamen-watercolor, white-bg  The built-
        in Mapbox styles are: basic, streets, outdoors, light, dark,
        satellite, satellite-streets  Mapbox style URLs are of the
        form: mapbox://mapbox.mapbox-<name>-<version>

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
        Sets the zoom level of the map (mapbox.zoom).

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
        accesstoken
            Sets the mapbox access token to be used for this mapbox
            map. Alternatively, the mapbox access token can be set
            in the configuration options under `mapboxAccessToken`.
            Note that accessToken are only required when `style`
            (e.g with values : basic, streets, outdoors, light,
            dark, satellite, satellite-streets ) and/or a layout
            layer references the Mapbox server.
        bearing
            Sets the bearing angle of the map in degrees counter-
            clockwise from North (mapbox.bearing).
        bounds
            :class:`plotly.graph_objects.layout.mapbox.Bounds`
            instance or dict with compatible properties
        center
            :class:`plotly.graph_objects.layout.mapbox.Center`
            instance or dict with compatible properties
        domain
            :class:`plotly.graph_objects.layout.mapbox.Domain`
            instance or dict with compatible properties
        layers
            A tuple of
            :class:`plotly.graph_objects.layout.mapbox.Layer`
            instances or dicts with compatible properties
        layerdefaults
            When used in a template (as
            layout.template.layout.mapbox.layerdefaults), sets the
            default property values to use for elements of
            layout.mapbox.layers
        pitch
            Sets the pitch angle of the map (in degrees, where 0
            means perpendicular to the surface of the map)
            (mapbox.pitch).
        style
            Defines the map layers that are rendered by default
            below the trace layers defined in `data`, which are
            themselves by default rendered below the layers defined
            in `layout.mapbox.layers`.  These layers can be defined
            either explicitly as a Mapbox Style object which can
            contain multiple layer definitions that load data from
            any public or private Tile Map Service (TMS or XYZ) or
            Web Map Service (WMS) or implicitly by using one of the
            built-in style objects which use WMSes which do not
            require any access tokens, or by using a default Mapbox
            style or custom Mapbox style URL, both of which require
            a Mapbox access token  Note that Mapbox access token
            can be set in the `accesstoken` attribute or in the
            `mapboxAccessToken` config option.  Mapbox Style
            objects are of the form described in the Mapbox GL JS
            documentation available at
            https://docs.mapbox.com/mapbox-gl-js/style-spec  The
            built-in plotly.js styles objects are: carto-
            darkmatter, carto-positron, open-street-map, stamen-
            terrain, stamen-toner, stamen-watercolor, white-bg  The
            built-in Mapbox styles are: basic, streets, outdoors,
            light, dark, satellite, satellite-streets  Mapbox style
            URLs are of the form:
            mapbox://mapbox.mapbox-<name>-<version>
        uirevision
            Controls persistence of user-driven changes in the
            view: `center`, `zoom`, `bearing`, `pitch`. Defaults to
            `layout.uirevision`.
        zoom
            Sets the zoom level of the map (mapbox.zoom).
        """

    def __init__(
        self,
        arg=None,
        accesstoken=None,
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
        Construct a new Mapbox object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Mapbox`
        accesstoken
            Sets the mapbox access token to be used for this mapbox
            map. Alternatively, the mapbox access token can be set
            in the configuration options under `mapboxAccessToken`.
            Note that accessToken are only required when `style`
            (e.g with values : basic, streets, outdoors, light,
            dark, satellite, satellite-streets ) and/or a layout
            layer references the Mapbox server.
        bearing
            Sets the bearing angle of the map in degrees counter-
            clockwise from North (mapbox.bearing).
        bounds
            :class:`plotly.graph_objects.layout.mapbox.Bounds`
            instance or dict with compatible properties
        center
            :class:`plotly.graph_objects.layout.mapbox.Center`
            instance or dict with compatible properties
        domain
            :class:`plotly.graph_objects.layout.mapbox.Domain`
            instance or dict with compatible properties
        layers
            A tuple of
            :class:`plotly.graph_objects.layout.mapbox.Layer`
            instances or dicts with compatible properties
        layerdefaults
            When used in a template (as
            layout.template.layout.mapbox.layerdefaults), sets the
            default property values to use for elements of
            layout.mapbox.layers
        pitch
            Sets the pitch angle of the map (in degrees, where 0
            means perpendicular to the surface of the map)
            (mapbox.pitch).
        style
            Defines the map layers that are rendered by default
            below the trace layers defined in `data`, which are
            themselves by default rendered below the layers defined
            in `layout.mapbox.layers`.  These layers can be defined
            either explicitly as a Mapbox Style object which can
            contain multiple layer definitions that load data from
            any public or private Tile Map Service (TMS or XYZ) or
            Web Map Service (WMS) or implicitly by using one of the
            built-in style objects which use WMSes which do not
            require any access tokens, or by using a default Mapbox
            style or custom Mapbox style URL, both of which require
            a Mapbox access token  Note that Mapbox access token
            can be set in the `accesstoken` attribute or in the
            `mapboxAccessToken` config option.  Mapbox Style
            objects are of the form described in the Mapbox GL JS
            documentation available at
            https://docs.mapbox.com/mapbox-gl-js/style-spec  The
            built-in plotly.js styles objects are: carto-
            darkmatter, carto-positron, open-street-map, stamen-
            terrain, stamen-toner, stamen-watercolor, white-bg  The
            built-in Mapbox styles are: basic, streets, outdoors,
            light, dark, satellite, satellite-streets  Mapbox style
            URLs are of the form:
            mapbox://mapbox.mapbox-<name>-<version>
        uirevision
            Controls persistence of user-driven changes in the
            view: `center`, `zoom`, `bearing`, `pitch`. Defaults to
            `layout.uirevision`.
        zoom
            Sets the zoom level of the map (mapbox.zoom).

        Returns
        -------
        Mapbox
        """
        super(Mapbox, self).__init__("mapbox")

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
The first argument to the plotly.graph_objs.layout.Mapbox
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Mapbox`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("accesstoken", None)
        _v = accesstoken if accesstoken is not None else _v
        if _v is not None:
            self["accesstoken"] = _v
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
