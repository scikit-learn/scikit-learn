import _plotly_utils.basevalidators


class LabelsuffixValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="labelsuffix", parent_name="carpet.aaxis", **kwargs):
        super(LabelsuffixValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
