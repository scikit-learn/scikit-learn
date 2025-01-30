import _plotly_utils.basevalidators


class LeafValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="leaf", parent_name="sunburst", **kwargs):
        super(LeafValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Leaf"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            opacity
                Sets the opacity of the leaves. With colorscale
                it is defaulted to 1; otherwise it is defaulted
                to 0.7
""",
            ),
            **kwargs,
        )
