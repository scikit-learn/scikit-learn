import _plotly_utils.basevalidators


class ZerolinewidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="zerolinewidth", parent_name="layout.yaxis", **kwargs
    ):
        super(ZerolinewidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            **kwargs,
        )
