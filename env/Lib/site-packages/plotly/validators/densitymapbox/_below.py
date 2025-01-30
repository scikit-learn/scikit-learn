import _plotly_utils.basevalidators


class BelowValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="below", parent_name="densitymapbox", **kwargs):
        super(BelowValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
