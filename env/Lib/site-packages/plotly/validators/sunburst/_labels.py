import _plotly_utils.basevalidators


class LabelsValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="labels", parent_name="sunburst", **kwargs):
        super(LabelsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
