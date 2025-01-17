import _plotly_utils.basevalidators


class YbingroupValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="ybingroup", parent_name="histogram2d", **kwargs):
        super(YbingroupValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
