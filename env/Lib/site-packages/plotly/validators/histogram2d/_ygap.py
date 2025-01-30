import _plotly_utils.basevalidators


class YgapValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="ygap", parent_name="histogram2d", **kwargs):
        super(YgapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
