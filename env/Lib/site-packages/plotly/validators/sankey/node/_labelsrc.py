import _plotly_utils.basevalidators


class LabelsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="labelsrc", parent_name="sankey.node", **kwargs):
        super(LabelsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
