import _plotly_utils.basevalidators


class RootValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="root", parent_name="treemap", **kwargs):
        super(RootValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Root"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                sets the color of the root node for a
                sunburst/treemap/icicle trace. this has no
                effect when a colorscale is used to set the
                markers.
""",
            ),
            **kwargs,
        )
