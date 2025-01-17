import _plotly_utils.basevalidators


class LegendgrouptitleValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="legendgrouptitle", parent_name="scattergl", **kwargs
    ):
        super(LegendgrouptitleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Legendgrouptitle"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Sets this legend group's title font.
            text
                Sets the title of the legend group.
""",
            ),
            **kwargs,
        )
