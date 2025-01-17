import _plotly_utils.basevalidators


class StepmodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self,
        plotly_name="stepmode",
        parent_name="layout.xaxis.rangeselector.button",
        **kwargs,
    ):
        super(StepmodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["backward", "todate"]),
            **kwargs,
        )
