import _plotly_utils.basevalidators


class XperiodalignmentValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="xperiodalignment", parent_name="heatmap", **kwargs):
        super(XperiodalignmentValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            implied_edits=kwargs.pop("implied_edits", {"xtype": "scaled"}),
            values=kwargs.pop("values", ["start", "middle", "end"]),
            **kwargs,
        )
