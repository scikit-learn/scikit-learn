import _plotly_utils.basevalidators


class ProjectionValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="projection", parent_name="scatter3d", **kwargs):
        super(ProjectionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Projection"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            x
                :class:`plotly.graph_objects.scatter3d.projecti
                on.X` instance or dict with compatible
                properties
            y
                :class:`plotly.graph_objects.scatter3d.projecti
                on.Y` instance or dict with compatible
                properties
            z
                :class:`plotly.graph_objects.scatter3d.projecti
                on.Z` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
