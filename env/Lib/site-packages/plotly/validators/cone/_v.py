import _plotly_utils.basevalidators


class VValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="v", parent_name="cone", **kwargs):
        super(VValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
