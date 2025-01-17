import _plotly_utils.basevalidators


class HighValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="high", parent_name="candlestick", **kwargs):
        super(HighValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
