import _plotly_utils.basevalidators


class AddValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="add", parent_name="layout.modebar", **kwargs):
        super(AddValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "modebar"),
            **kwargs,
        )
