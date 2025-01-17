import _plotly_utils.basevalidators


class CountValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self,
        plotly_name="count",
        parent_name="layout.xaxis.rangeselector.button",
        **kwargs,
    ):
        super(CountValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
