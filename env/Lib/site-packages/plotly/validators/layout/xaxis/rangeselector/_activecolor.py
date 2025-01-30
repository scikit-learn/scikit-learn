import _plotly_utils.basevalidators


class ActivecolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(
        self,
        plotly_name="activecolor",
        parent_name="layout.xaxis.rangeselector",
        **kwargs,
    ):
        super(ActivecolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
