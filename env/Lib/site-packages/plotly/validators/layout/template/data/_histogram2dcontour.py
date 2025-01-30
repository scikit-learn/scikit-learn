import _plotly_utils.basevalidators


class Histogram2DcontourValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(
        self,
        plotly_name="histogram2dcontour",
        parent_name="layout.template.data",
        **kwargs,
    ):
        super(Histogram2DcontourValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Histogram2dContour"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
