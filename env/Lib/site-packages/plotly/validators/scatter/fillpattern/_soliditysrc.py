import _plotly_utils.basevalidators


class SoliditysrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="soliditysrc", parent_name="scatter.fillpattern", **kwargs
    ):
        super(SoliditysrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
