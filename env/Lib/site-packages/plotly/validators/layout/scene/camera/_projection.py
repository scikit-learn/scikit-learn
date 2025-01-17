import _plotly_utils.basevalidators


class ProjectionValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="projection", parent_name="layout.scene.camera", **kwargs
    ):
        super(ProjectionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Projection"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            type
                Sets the projection type. The projection type
                could be either "perspective" or
                "orthographic". The default is "perspective".
""",
            ),
            **kwargs,
        )
