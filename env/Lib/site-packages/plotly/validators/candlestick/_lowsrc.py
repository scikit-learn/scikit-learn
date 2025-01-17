import _plotly_utils.basevalidators


class LowsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="lowsrc", parent_name="candlestick", **kwargs):
        super(LowsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
