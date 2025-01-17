import _plotly_utils.basevalidators


class ProjectionValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="projection", parent_name="layout.geo", **kwargs):
        super(ProjectionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Projection"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            distance
                For satellite projection type only. Sets the
                distance from the center of the sphere to the
                point of view as a proportion of the sphere’s
                radius.
            parallels
                For conic projection types only. Sets the
                parallels (tangent, secant) where the cone
                intersects the sphere.
            rotation
                :class:`plotly.graph_objects.layout.geo.project
                ion.Rotation` instance or dict with compatible
                properties
            scale
                Zooms in or out on the map view. A scale of 1
                corresponds to the largest zoom level that fits
                the map's lon and lat ranges.
            tilt
                For satellite projection type only. Sets the
                tilt angle of perspective projection.
            type
                Sets the projection type.
""",
            ),
            **kwargs,
        )
