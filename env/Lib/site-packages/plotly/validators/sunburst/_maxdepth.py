import _plotly_utils.basevalidators


class MaxdepthValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="maxdepth", parent_name="sunburst", **kwargs):
        super(MaxdepthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
