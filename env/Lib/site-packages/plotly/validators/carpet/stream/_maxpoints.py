import _plotly_utils.basevalidators


class MaxpointsValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="maxpoints", parent_name="carpet.stream", **kwargs):
        super(MaxpointsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 10000),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
