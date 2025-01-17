import _plotly_utils.basevalidators


class XValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="x", parent_name="scatter3d.projection", **kwargs):
        super(XValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "X"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            opacity
                Sets the projection color.
            scale
                Sets the scale factor determining the size of
                the projection marker points.
            show
                Sets whether or not projections are shown along
                the x axis.
""",
            ),
            **kwargs,
        )
