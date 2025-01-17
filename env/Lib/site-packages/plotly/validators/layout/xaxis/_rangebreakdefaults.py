import _plotly_utils.basevalidators


class RangebreakdefaultsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="rangebreakdefaults", parent_name="layout.xaxis", **kwargs
    ):
        super(RangebreakdefaultsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Rangebreak"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
