import _plotly_utils.basevalidators


class EditableValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="editable", parent_name="layout.shape", **kwargs):
        super(EditableValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+arraydraw"),
            **kwargs,
        )
