import _plotly_utils.basevalidators


class TicklabelindexValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self, plotly_name="ticklabelindex", parent_name="layout.xaxis", **kwargs
    ):
        super(TicklabelindexValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
