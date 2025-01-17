import _plotly_utils.basevalidators


class TicklabelstandoffValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self, plotly_name="ticklabelstandoff", parent_name="layout.xaxis", **kwargs
    ):
        super(TicklabelstandoffValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            **kwargs,
        )
