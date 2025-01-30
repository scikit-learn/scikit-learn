import _plotly_utils.basevalidators


class CountValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="count", parent_name="volume.surface", **kwargs):
        super(CountValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
