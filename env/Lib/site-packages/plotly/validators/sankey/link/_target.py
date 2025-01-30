import _plotly_utils.basevalidators


class TargetValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="target", parent_name="sankey.link", **kwargs):
        super(TargetValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
