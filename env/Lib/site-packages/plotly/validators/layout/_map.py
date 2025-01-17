import _plotly_utils.basevalidators


class MapValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="map", parent_name="layout", **kwargs):
        super(MapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Map"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            bearing
                Sets the bearing angle of the map in degrees
                counter-clockwise from North (map.bearing).
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
                layout.template.layout.map.layerdefaults), sets
                the default property values to use for elements
                of layout.map.layers
            pitch
                Sets the pitch angle of the map (in degrees,
                where 0 means perpendicular to the surface of
                the map) (map.pitch).
            style
                Defines the map layers that are rendered by
                default below the trace layers defined in
                `data`, which are themselves by default
                rendered below the layers defined in
                `layout.map.layers`.  These layers can be
                defined either explicitly as a Map Style object
                which can contain multiple layer definitions
                that load data from any public or private Tile
                Map Service (TMS or XYZ) or Web Map Service
                (WMS) or implicitly by using one of the built-
                in style objects which use WMSes or by using a
                custom style URL  Map Style objects are of the
                form described in the MapLibre GL JS
                documentation available at
                https://maplibre.org/maplibre-style-spec/  The
                built-in plotly.js styles objects are: basic,
                carto-darkmatter, carto-darkmatter-nolabels,
                carto-positron, carto-positron-nolabels, carto-
                voyager, carto-voyager-nolabels, dark, light,
                open-street-map, outdoors, satellite,
                satellite-streets, streets, white-bg.
            uirevision
                Controls persistence of user-driven changes in
                the view: `center`, `zoom`, `bearing`, `pitch`.
                Defaults to `layout.uirevision`.
            zoom
                Sets the zoom level of the map (map.zoom).
""",
            ),
            **kwargs,
        )
