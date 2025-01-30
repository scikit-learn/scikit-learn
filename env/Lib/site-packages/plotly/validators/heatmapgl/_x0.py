import _plotly_utils.basevalidators


class X0Validator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="x0", parent_name="heatmapgl", **kwargs):
        super(X0Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            implied_edits=kwargs.pop("implied_edits", {"xtype": "scaled"}),
            **kwargs,
        )
