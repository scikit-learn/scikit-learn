import _plotly_utils.basevalidators


class CenterValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="center", parent_name="layout.geo", **kwargs):
        super(CenterValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Center"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            lat
                Sets the latitude of the map's center. For all
                projection types, the map's latitude center
                lies at the middle of the latitude range by
                default.
            lon
                Sets the longitude of the map's center. By
                default, the map's longitude center lies at the
                middle of the longitude range for scoped
                projection and above `projection.rotation.lon`
                otherwise.
""",
            ),
            **kwargs,
        )
