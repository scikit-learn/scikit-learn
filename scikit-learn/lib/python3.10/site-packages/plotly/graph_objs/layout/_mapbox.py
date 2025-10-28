#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Mapbox(_BaseLayoutHierarchyType):
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

    @property
    def bounds(self):
        """
        The 'bounds' property is an instance of Bounds
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.mapbox.Bounds`
          - A dict of string/value properties that will be passed
            to the Bounds constructor

        Returns
        -------
        plotly.graph_objs.layout.mapbox.Bounds
        """
        return self["bounds"]

    @bounds.setter
    def bounds(self, val):
        self["bounds"] = val

    @property
    def center(self):
        """
        The 'center' property is an instance of Center
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.mapbox.Center`
          - A dict of string/value properties that will be passed
            to the Center constructor

        Returns
        -------
        plotly.graph_objs.layout.mapbox.Center
        """
        return self["center"]

    @center.setter
    def center(self, val):
        self["center"] = val

    @property
    def domain(self):
        """
        The 'domain' property is an instance of Domain
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.mapbox.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

        Returns
        -------
        plotly.graph_objs.layout.mapbox.Domain
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    @property
    def layers(self):
        """
        The 'layers' property is a tuple of instances of
        Layer that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.mapbox.Layer
          - A list or tuple of dicts of string/value properties that
            will be passed to the Layer constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.mapbox.Layer]
        """
        return self["layers"]

    @layers.setter
    def layers(self, val):
        self["layers"] = val

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

        Returns
        -------
        plotly.graph_objs.layout.mapbox.Layer
        """
        return self["layerdefaults"]

    @layerdefaults.setter
    def layerdefaults(self, val):
        self["layerdefaults"] = val

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
        super().__init__("mapbox")
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
The first argument to the plotly.graph_objs.layout.Mapbox
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Mapbox`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("accesstoken", arg, accesstoken)
        self._set_property("bearing", arg, bearing)
        self._set_property("bounds", arg, bounds)
        self._set_property("center", arg, center)
        self._set_property("domain", arg, domain)
        self._set_property("layers", arg, layers)
        self._set_property("layerdefaults", arg, layerdefaults)
        self._set_property("pitch", arg, pitch)
        self._set_property("style", arg, style)
        self._set_property("uirevision", arg, uirevision)
        self._set_property("zoom", arg, zoom)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
