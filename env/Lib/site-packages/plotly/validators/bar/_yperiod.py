import _plotly_utils.basevalidators


class YperiodValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="yperiod", parent_name="bar", **kwargs):
        super(YperiodValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
