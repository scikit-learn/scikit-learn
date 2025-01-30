import _plotly_utils.basevalidators


class DiagonalValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="diagonal", parent_name="splom", **kwargs):
        super(DiagonalValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Diagonal"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            visible
                Determines whether or not subplots on the
                diagonal are displayed.
""",
            ),
            **kwargs,
        )
