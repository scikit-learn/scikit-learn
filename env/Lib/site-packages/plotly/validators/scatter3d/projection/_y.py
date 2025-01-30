import _plotly_utils.basevalidators


class YValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="y", parent_name="scatter3d.projection", **kwargs):
        super(YValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Y"),
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
                the y axis.
""",
            ),
            **kwargs,
        )
