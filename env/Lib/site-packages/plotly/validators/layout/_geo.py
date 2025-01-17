import _plotly_utils.basevalidators


class GeoValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="geo", parent_name="layout", **kwargs):
        super(GeoValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Geo"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            bgcolor
                Set the background color of the map
            center
                :class:`plotly.graph_objects.layout.geo.Center`
                instance or dict with compatible properties
            coastlinecolor
                Sets the coastline color.
            coastlinewidth
                Sets the coastline stroke width (in px).
            countrycolor
                Sets line color of the country boundaries.
            countrywidth
                Sets line width (in px) of the country
                boundaries.
            domain
                :class:`plotly.graph_objects.layout.geo.Domain`
                instance or dict with compatible properties
            fitbounds
                Determines if this subplot's view settings are
                auto-computed to fit trace data. On scoped
                maps, setting `fitbounds` leads to `center.lon`
                and `center.lat` getting auto-filled. On maps
                with a non-clipped projection, setting
                `fitbounds` leads to `center.lon`,
                `center.lat`, and `projection.rotation.lon`
                getting auto-filled. On maps with a clipped
                projection, setting `fitbounds` leads to
                `center.lon`, `center.lat`,
                `projection.rotation.lon`,
                `projection.rotation.lat`, `lonaxis.range` and
                `lonaxis.range` getting auto-filled. If
                "locations", only the trace's visible locations
                are considered in the `fitbounds` computations.
                If "geojson", the entire trace input `geojson`
                (if provided) is considered in the `fitbounds`
                computations, Defaults to False.
            framecolor
                Sets the color the frame.
            framewidth
                Sets the stroke width (in px) of the frame.
            lakecolor
                Sets the color of the lakes.
            landcolor
                Sets the land mass color.
            lataxis
                :class:`plotly.graph_objects.layout.geo.Lataxis
                ` instance or dict with compatible properties
            lonaxis
                :class:`plotly.graph_objects.layout.geo.Lonaxis
                ` instance or dict with compatible properties
            oceancolor
                Sets the ocean color
            projection
                :class:`plotly.graph_objects.layout.geo.Project
                ion` instance or dict with compatible
                properties
            resolution
                Sets the resolution of the base layers. The
                values have units of km/mm e.g. 110 corresponds
                to a scale ratio of 1:110,000,000.
            rivercolor
                Sets color of the rivers.
            riverwidth
                Sets the stroke width (in px) of the rivers.
            scope
                Set the scope of the map.
            showcoastlines
                Sets whether or not the coastlines are drawn.
            showcountries
                Sets whether or not country boundaries are
                drawn.
            showframe
                Sets whether or not a frame is drawn around the
                map.
            showlakes
                Sets whether or not lakes are drawn.
            showland
                Sets whether or not land masses are filled in
                color.
            showocean
                Sets whether or not oceans are filled in color.
            showrivers
                Sets whether or not rivers are drawn.
            showsubunits
                Sets whether or not boundaries of subunits
                within countries (e.g. states, provinces) are
                drawn.
            subunitcolor
                Sets the color of the subunits boundaries.
            subunitwidth
                Sets the stroke width (in px) of the subunits
                boundaries.
            uirevision
                Controls persistence of user-driven changes in
                the view (projection and center). Defaults to
                `layout.uirevision`.
            visible
                Sets the default visibility of the base layers.
""",
            ),
            **kwargs,
        )
