import _plotly_utils.basevalidators


class DecreasingValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="decreasing", parent_name="indicator.delta", **kwargs
    ):
        super(DecreasingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Decreasing"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the color for increasing value.
            symbol
                Sets the symbol to display for increasing value
""",
            ),
            **kwargs,
        )
